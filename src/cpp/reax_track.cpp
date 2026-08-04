#include "reax_track.h"

#include <climits>
#include <cmath>
#include <filesystem>
#include <functional>
#include <queue>
#include <sstream>

#include "fmt/format.h"
#include "argparser.h"
#include "string_tools.h"

// External global variables (defined in argparser.cpp)
extern bool FLAG_TRACK_REACTIONS;
extern int STABLE_TIME_FRAMES;
extern float TIMESTEP_FS;
extern int SAMPLING_FREQ;

static unsigned int molecule_identity_id(const std::string& formula) {
    return reax_string_hash(formula);
}

static std::string json_escape(const std::string& value) {
    std::string out;
    for (char c : value) {
        if (c == '"' || c == '\\') {
            out += '\\';
            out += c;
        } else if (c == '\n') {
            out += "\\n";
        } else if (c == '\r') {
            out += "\\r";
        } else if (c == '\t') {
            out += "\\t";
        } else {
            out += c;
        }
    }
    return out;
}

static std::string bond_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return std::to_string(a) + "-" + std::to_string(b);
}

static float coord_distance(const std::vector<float>& a, const std::vector<float>& b) {
    float dx = a[0] - b[0];
    float dy = a[1] - b[1];
    float dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static std::vector<float> minimum_image_delta(const std::vector<float>& from,
                                              const std::vector<float>& to,
                                              bool has_boundaries,
                                              const std::vector<float>& axis_lengths) {
    std::vector<float> delta(3, 0.0f);
    for (size_t i = 0; i < 3; ++i) {
        delta[i] = to[i] - from[i];
        if (has_boundaries && i < axis_lengths.size() && axis_lengths[i] > 0.0f) {
            delta[i] = delta[i] - axis_lengths[i] * std::floor(delta[i] / axis_lengths[i] + 0.5f);
        }
    }
    return delta;
}

static float coord_distance_snapshot(const AtomSnapshot& a, const AtomSnapshot& b) {
    float dx = a.x - b.x;
    float dy = a.y - b.y;
    float dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static bool cancel_common_species(std::vector<std::string>& reactants,
                                    std::vector<std::string>& products) {
    std::map<std::string, int> left_counts, right_counts;
    for (const auto& s : reactants) left_counts[s]++;
    for (const auto& s : products) right_counts[s]++;

    for (auto& [species, left_n] : left_counts) {
        auto it = right_counts.find(species);
        if (it != right_counts.end()) {
            int cancel_n = std::min(left_n, it->second);
            left_n -= cancel_n;
            it->second -= cancel_n;
        }
    }

    reactants.clear();
    products.clear();
    for (const auto& [s, n] : left_counts) {
        for (int i = 0; i < n; i++) reactants.push_back(s);
    }
    for (const auto& [s, n] : right_counts) {
        for (int i = 0; i < n; i++) products.push_back(s);
    }

    return !reactants.empty() || !products.empty();
}

static bool check_formula_atom_balance(const std::vector<std::string>& reactants,
                                        const std::vector<std::string>& products) {
    std::map<std::string, int> total;
    for (const auto& formula : reactants) {
        for (const auto& [elem, count] : parse_formula(formula)) {
            total[elem] += count;
        }
    }
    for (const auto& formula : products) {
        for (const auto& [elem, count] : parse_formula(formula)) {
            total[elem] -= count;
        }
    }
    for (const auto& [elem, count] : total) {
        if (count != 0) return false;
    }
    return true;
}

// Static member initialization
int TrackedMolecule::next_id = 0;

//=============================================================================
// TrackedMolecule Implementation
//=============================================================================

TrackedMolecule::TrackedMolecule(int frame, const Molecule* mol, bool has_boundaries, const std::vector<float>& axis_lengths)
    : id(next_id++),
      start_frame(frame),
      end_frame(INT_MAX),
      formula(mol->formula),
      hash(mol->hash),
      is_stable(false),
      is_processed(false) {
    
    // Copy atom IDs
    for (int atom_id : mol->atom_ids) {
        atom_ids.insert(atom_id);
    }
    std::map<int, const Atom*> atom_by_id;
    for (const auto* atom : mol->mol_atoms) {
        atom_by_id[atom->id] = atom;
    }

    std::map<int, std::vector<int>> graph;
    for (const auto* bond : mol->mol_bonds) {
        int a = bond->atom_i->id;
        int b = bond->atom_j->id;
        graph[a].push_back(b);
        graph[b].push_back(a);
    }

    std::map<int, std::vector<float>> unwrapped;
    if (!atom_by_id.empty()) {
        int root = atom_by_id.begin()->first;
        unwrapped[root] = atom_by_id[root]->coord;
        std::queue<int> queue;
        queue.push(root);
        while (!queue.empty()) {
            int current = queue.front();
            queue.pop();
            const auto& current_raw = atom_by_id[current]->coord;
            const auto& current_unwrapped = unwrapped[current];
            for (int neighbor : graph[current]) {
                if (unwrapped.count(neighbor)) continue;
                auto delta = minimum_image_delta(current_raw, atom_by_id[neighbor]->coord, has_boundaries, axis_lengths);
                unwrapped[neighbor] = {
                    current_unwrapped[0] + delta[0],
                    current_unwrapped[1] + delta[1],
                    current_unwrapped[2] + delta[2],
                };
                queue.push(neighbor);
            }
        }
        for (const auto& [atom_id, atom] : atom_by_id) {
            if (!unwrapped.count(atom_id)) {
                unwrapped[atom_id] = atom->coord;
            }
        }
    }

    atoms.reserve(mol->mol_atoms.size());
    for (const auto* atom : mol->mol_atoms) {
        const auto& coord = unwrapped[atom->id];
        atoms.push_back({atom->id, atom->type_name, coord[0], coord[1], coord[2]});
    }
    std::sort(atoms.begin(), atoms.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
    std::map<int, AtomSnapshot> snapshot_by_id;
    for (const auto& atom : atoms) {
        snapshot_by_id[atom.id] = atom;
    }

    bonds.reserve(mol->mol_bonds.size());
    for (const auto* bond : mol->mol_bonds) {
        int a = bond->atom_i->id;
        int b = bond->atom_j->id;
        if (a > b) std::swap(a, b);
        bonds.push_back({a, b, bond->order, coord_distance_snapshot(snapshot_by_id[a], snapshot_by_id[b])});
    }
    std::sort(bonds.begin(), bonds.end(), [](const auto& a, const auto& b) {
        if (a.atom_i != b.atom_i) return a.atom_i < b.atom_i;
        return a.atom_j < b.atom_j;
    });
}

std::string TrackedMolecule::to_string() const {
    return fmt::format("{}(id={},frames={}-{})", formula, id, start_frame, 
                       is_active() ? "inf" : std::to_string(end_frame));
}

//=============================================================================
// ReactionEvent Implementation
//=============================================================================

ReactionEvent::ReactionEvent(int eid, int fid, float time)
    : event_id(eid),
      frame_id(fid),
      time_fs(time) {}

std::string ReactionEvent::get_reactants_string() const {
    std::string result;
    for (size_t i = 0; i < reactants.size(); ++i) {
        if (i > 0) result += "+";
        result += reactants[i]->formula;
    }
    return result;
}

std::string ReactionEvent::get_products_string() const {
    std::string result;
    for (size_t i = 0; i < products.size(); ++i) {
        if (i > 0) result += "+";
        result += products[i]->formula;
    }
    return result;
}

bool ReactionEvent::check_atom_conservation() const {
    std::unordered_set<int> reactant_atoms;
    std::unordered_set<int> product_atoms;
    
    for (const auto* mol : reactants) {
        reactant_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    }
    for (const auto* mol : products) {
        product_atoms.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    }
    
    return reactant_atoms == product_atoms;
}

static std::string atom_ids_string(const std::vector<TrackedMolecule*>& molecules) {
    std::vector<int> atom_ids;
    for (const auto* mol : molecules) {
        atom_ids.insert(atom_ids.end(), mol->atom_ids.begin(), mol->atom_ids.end());
    }
    std::sort(atom_ids.begin(), atom_ids.end());
    atom_ids.erase(std::unique(atom_ids.begin(), atom_ids.end()), atom_ids.end());

    std::string result;
    for (size_t i = 0; i < atom_ids.size(); ++i) {
        if (i > 0) result += "+";
        result += std::to_string(atom_ids[i]);
    }
    return result;
}

//=============================================================================
// ReactionTracker Implementation
//=============================================================================

ReactionTracker::ReactionTracker(int stable_frames, float ts, int freq)
    : stable_time_frames(stable_frames),
      timestep_fs(ts),
      sampling_frequency(freq),
      next_event_id(0),
      last_processed_frame(0) {}

ReactionTracker::~ReactionTracker() {
    // Clean up all tracked molecules
    for (auto* mol : all_molecules) {
        delete mol;
    }
    all_molecules.clear();
}

TrackedMolecule* ReactionTracker::find_or_create_tracked(const Molecule* mol, int frame_id, bool has_boundaries, const std::vector<float>& axis_lengths) {
    // First, check if any active molecule has the same atom set and structure.
    for (auto* tracked : active_molecules) {
        if (tracked->atom_ids == mol->atom_ids && tracked->hash == mol->hash) {
            return tracked;
        }
    }
    
    // Create new tracked molecule
    TrackedMolecule* tracked = new TrackedMolecule(frame_id, mol, has_boundaries, axis_lengths);
    all_molecules.push_back(tracked);
    active_molecules.insert(tracked);
    
    return tracked;
}

bool ReactionTracker::should_mark_stable(TrackedMolecule* mol, int current_frame) {
    if (mol->is_stable) return false;
    if (mol->is_active()) {
        return (current_frame - mol->start_frame) >= stable_time_frames;
    }
    return (mol->end_frame - mol->start_frame) >= stable_time_frames;
}

void ReactionTracker::create_reaction_event(int frame_id,
                                           const std::vector<TrackedMolecule*>& reactants,
                                           const std::vector<TrackedMolecule*>& products) {
    // Calculate time in fs
    float time = frame_id * timestep_fs * sampling_frequency;
    
    ReactionEvent event(next_event_id++, frame_id, time);
    event.reactants = reactants;
    event.products = products;
    
    // Verify atom conservation
    if (!event.check_atom_conservation()) {
        fmt::print("Warning: Reaction event {} does not conserve atoms!\n", event.event_id);
    }
    
    events.push_back(std::move(event));
}

void ReactionTracker::detect_reactions_from_changes(
        int frame_id,
        const std::unordered_map<int, Molecule*>& current_atom_to_mol,
        const std::unordered_set<int>& changed_atoms,
        bool has_boundaries,
        const std::vector<float>& axis_lengths) {
    
    std::unordered_set<int> processed_atoms;
    
    for (int start_atom : changed_atoms) {
        if (processed_atoms.count(start_atom)) continue;
        
        // Symmetric difference algorithm to find complete reaction cluster
        std::unordered_set<int> sym_diff = {start_atom};
        std::vector<TrackedMolecule*> reactants;
        std::vector<TrackedMolecule*> products;
        
        std::unordered_set<TrackedMolecule*> used_reactants;
        std::unordered_set<TrackedMolecule*> used_products;
        
        while (!sym_diff.empty()) {
            int atom_id = *sym_diff.begin();
            sym_diff.erase(sym_diff.begin());
            
            if (processed_atoms.count(atom_id)) continue;
            processed_atoms.insert(atom_id);
            
            // Check old molecule (from previous frame)
            auto old_it = atom_to_molecule.find(atom_id);
            if (old_it != atom_to_molecule.end()) {
                TrackedMolecule* old_tracked = old_it->second;
                if (!used_reactants.count(old_tracked)) {
                    reactants.push_back(old_tracked);
                    used_reactants.insert(old_tracked);
                    
                    // Add all atoms from this molecule to sym_diff
                    for (int a : old_tracked->atom_ids) {
                        if (!processed_atoms.count(a)) {
                            sym_diff.insert(a);
                        }
                    }
                }
            }
            
            // Check new molecule (in current frame)
            auto new_it = current_atom_to_mol.find(atom_id);
            if (new_it != current_atom_to_mol.end()) {
                Molecule* new_mol = new_it->second;
                TrackedMolecule* new_tracked = find_or_create_tracked(new_mol, frame_id, has_boundaries, axis_lengths);
                
                if (!used_products.count(new_tracked)) {
                    products.push_back(new_tracked);
                    used_products.insert(new_tracked);
                    
                    // Add all atoms from this molecule to sym_diff
                    for (int a : new_tracked->atom_ids) {
                        if (!processed_atoms.count(a)) {
                            sym_diff.insert(a);
                        }
                    }
                }
            }
        }
        
        // Only record if we have both reactants and products
        if (!reactants.empty() && !products.empty()) {
            // Set up predecessor/successor relationships
            for (auto* r : reactants) {
                for (auto* p : products) {
                    r->successors.push_back(p);
                    p->predecessors.push_back(r);
                }
            }
            
            // Mark reactants as ended
            for (auto* r : reactants) {
                if (r->is_active()) {
                    r->end_frame = frame_id;
                }
            }
            
            // Create reaction event
            create_reaction_event(frame_id, reactants, products);
        }
    }
}

void ReactionTracker::process_frame(int frame_id, const std::vector<Molecule*>& molecules, bool has_boundaries, const std::vector<float>& axis_lengths) {
    if (!FLAG_TRACK_REACTIONS) return;
    
    // Build atom -> molecule mapping for current frame
    std::unordered_map<int, Molecule*> current_atom_to_mol;
    for (Molecule* mol : molecules) {
        for (int atom_id : mol->atom_ids) {
            current_atom_to_mol[atom_id] = mol;
        }
    }
    
    // First frame: just initialize mappings
    if (frame_id == 1 || atom_to_molecule.empty()) {
        for (Molecule* mol : molecules) {
            TrackedMolecule* tracked = find_or_create_tracked(mol, frame_id, has_boundaries, axis_lengths);
            for (int atom_id : mol->atom_ids) {
                atom_to_molecule[atom_id] = tracked;
            }
        }
        last_processed_frame = frame_id;
        return;
    }
    
    // Find atoms that changed molecular affiliation
    std::unordered_set<int> changed_atoms;
    
    for (const auto& [atom_id, new_mol] : current_atom_to_mol) {
        auto old_it = atom_to_molecule.find(atom_id);
        if (old_it != atom_to_molecule.end()) {
            TrackedMolecule* old_tracked = old_it->second;
            // Check if molecule hash changed
            if (old_tracked->hash != new_mol->hash) {
                changed_atoms.insert(atom_id);
            }
        } else {
            // New atom appeared (shouldn't happen in standard MD)
            changed_atoms.insert(atom_id);
        }
    }
    
    // Detect reactions using symmetric difference algorithm
    if (!changed_atoms.empty()) {
        detect_reactions_from_changes(frame_id, current_atom_to_mol, changed_atoms, has_boundaries, axis_lengths);
    }
    
    // Update atom_to_molecule for next frame
    atom_to_molecule.clear();
    for (Molecule* mol : molecules) {
        TrackedMolecule* tracked = find_or_create_tracked(mol, frame_id, has_boundaries, axis_lengths);
        for (int atom_id : mol->atom_ids) {
            atom_to_molecule[atom_id] = tracked;
        }
    }
    
    // Mark molecules that are no longer present as ended
    std::unordered_set<TrackedMolecule*> current_mols;
    for (Molecule* mol : molecules) {
        current_mols.insert(find_or_create_tracked(mol, frame_id, has_boundaries, axis_lengths));
    }
    
    for (auto* tracked : active_molecules) {
        if (!current_mols.count(tracked) && tracked->is_active()) {
            tracked->end_frame = frame_id;
        }
    }
    
    // Update active molecules set
    active_molecules = std::move(current_mols);
    last_processed_frame = frame_id;
}

void ReactionTracker::finalize(int last_frame) {
    if (!FLAG_TRACK_REACTIONS) return;
    
    // Mark all remaining active molecules as ended
    for (auto* mol : active_molecules) {
        if (mol->is_active()) {
            mol->end_frame = last_frame;
        }
    }
}

void ReactionTracker::save_raw_events(const std::string& filepath) {
    if (events.empty()) return;

    FILE* fp = create_file(filepath);

    fmt::print(fp, "event_id,frame,n_reactants,n_products,reactant_hashes,product_hashes,reactant_formulas,product_formulas\n");

    for (const auto& event : events) {
        std::string reactant_hashes;
        std::string reactant_formulas;
        for (size_t i = 0; i < event.reactants.size(); ++i) {
            if (i > 0) reactant_hashes += "+";
            if (i > 0) reactant_formulas += "+";
            reactant_hashes += std::to_string(molecule_identity_id(event.reactants[i]->formula));
            reactant_formulas += event.reactants[i]->formula;
        }

        std::string product_hashes;
        std::string product_formulas;
        for (size_t i = 0; i < event.products.size(); ++i) {
            if (i > 0) product_hashes += "+";
            if (i > 0) product_formulas += "+";
            product_hashes += std::to_string(molecule_identity_id(event.products[i]->formula));
            product_formulas += event.products[i]->formula;
        }

        fmt::print(fp, "{},{},{},{},\"{}\",\"{}\",\"{}\",\"{}\"\n",
                   event.event_id,
                   event.frame_id,
                   event.reactants.size(),
                   event.products.size(),
                   reactant_hashes,
                   product_hashes,
                   reactant_formulas,
                   product_formulas);
    }

    fclose(fp);
    fmt::print("Saved {} raw reaction events to {}\n", events.size(), filepath);
}

static std::vector<std::string> sorted_hashes(const std::vector<TrackedMolecule*>& molecules) {
    std::vector<std::string> values;
    values.reserve(molecules.size());
    for (const auto* mol : molecules) {
        values.push_back(std::to_string(molecule_identity_id(mol->formula)));
    }
    std::sort(values.begin(), values.end());
    return values;
}

static std::string join_strings(const std::vector<std::string>& values, const std::string& sep = "+") {
    std::string out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) out += sep;
        out += values[i];
    }
    return out;
}

static std::set<int> atom_id_set(const std::vector<TrackedMolecule*>& molecules) {
    std::set<int> ids;
    for (const auto* mol : molecules) {
        ids.insert(mol->atom_ids.begin(), mol->atom_ids.end());
    }
    return ids;
}

static std::map<std::string, BondSnapshot> bond_map(const std::vector<TrackedMolecule*>& molecules) {
    std::map<std::string, BondSnapshot> out;
    for (const auto* mol : molecules) {
        for (const auto& bond : mol->bonds) {
            out[bond_key(bond.atom_i, bond.atom_j)] = bond;
        }
    }
    return out;
}

static void print_atoms_json(FILE* fp, const std::vector<TrackedMolecule*>& molecules, const char* indent) {
    bool first = true;
    for (const auto* mol : molecules) {
        for (const auto& atom : mol->atoms) {
            if (!first) fmt::print(fp, ",\n");
            first = false;
            fmt::print(fp,
                       "{}{{\"id\":{},\"element\":\"{}\",\"x\":{:.6f},\"y\":{:.6f},\"z\":{:.6f}}}",
                       indent,
                       atom.id,
                       json_escape(atom.element),
                       atom.x,
                       atom.y,
                       atom.z);
        }
    }
}

static void print_bonds_json(FILE* fp, const std::vector<TrackedMolecule*>& molecules, const char* indent) {
    bool first = true;
    for (const auto* mol : molecules) {
        for (const auto& bond : mol->bonds) {
            if (!first) fmt::print(fp, ",\n");
            first = false;
            fmt::print(fp,
                       "{}{{\"a\":{},\"b\":{},\"order\":{},\"length\":{:.6f}}}",
                       indent,
                       bond.atom_i,
                       bond.atom_j,
                       bond.order,
                       bond.length);
        }
    }
}

static void print_molecules_json(FILE* fp, const std::vector<TrackedMolecule*>& molecules, const char* indent) {
    for (size_t i = 0; i < molecules.size(); ++i) {
        const auto* mol = molecules[i];
        if (i > 0) fmt::print(fp, ",\n");
        std::vector<int> atom_ids(mol->atom_ids.begin(), mol->atom_ids.end());
        std::sort(atom_ids.begin(), atom_ids.end());
        fmt::print(fp,
                   "{}{{\"tracked_molecule_id\":{},\"molecule_hash\":\"{}\",\"formula\":\"{}\",\"atom_ids\":[",
                   indent,
                   mol->id,
                   molecule_identity_id(mol->formula),
                   json_escape(mol->formula));
        for (size_t j = 0; j < atom_ids.size(); ++j) {
            if (j > 0) fmt::print(fp, ",");
            fmt::print(fp, "{}", atom_ids[j]);
        }
        fmt::print(fp, "]}}");
    }
}

static void print_bond_change_ids_json(FILE* fp,
                                       const std::map<std::string, BondSnapshot>& lhs,
                                       const std::map<std::string, BondSnapshot>& rhs,
                                       const char* indent) {
    bool first = true;
    for (const auto& [key, bond] : lhs) {
        if (rhs.count(key)) continue;
        if (!first) fmt::print(fp, ",\n");
        first = false;
        fmt::print(fp, "{}\"{}\"", indent, key);
    }
}

void ReactionTracker::save_reaction_snapshots(const std::string& directory) {
    if (events.empty()) return;

    std::filesystem::create_directories(std::filesystem::path(OUTPUT_DIR) / directory);
    FILE* manifest = create_file(directory + "/snapshots_manifest.csv");
    fmt::print(manifest, "event_id,frame,signature,arity,path,atom_ids_conserved,reactant_atom_count,product_atom_count,broken_bond_count,created_bond_count\n");

    std::map<std::string, int> signature_counts;
    int written = 0;

    for (const auto& event : events) {
        if (event.reactants.empty() || event.products.empty()) continue;
        if (event.reactants.size() > 2 || event.products.size() > 2) continue;

        std::vector<std::string> reactant_hashes = sorted_hashes(event.reactants);
        std::vector<std::string> product_hashes = sorted_hashes(event.products);
        std::string signature = join_strings(reactant_hashes) + "->" + join_strings(product_hashes);
        if (signature_counts[signature] >= 3) continue;
        signature_counts[signature]++;

        auto reactant_atoms = atom_id_set(event.reactants);
        auto product_atoms = atom_id_set(event.products);
        auto reactant_bonds = bond_map(event.reactants);
        auto product_bonds = bond_map(event.products);
        bool atom_ids_conserved = reactant_atoms == product_atoms;

        int broken_count = 0;
        for (const auto& [key, _bond] : reactant_bonds) {
            if (!product_bonds.count(key)) broken_count++;
        }
        int created_count = 0;
        for (const auto& [key, _bond] : product_bonds) {
            if (!reactant_bonds.count(key)) created_count++;
        }

        std::ostringstream filename;
        filename << "event_" << event.event_id << ".rxtsnap.json";
        std::string rel_path = filename.str();
        std::string path = directory + "/" + rel_path;
        FILE* fp = create_file(path);

        fmt::print(fp, "{{\n");
        fmt::print(fp, "  \"schema\": \"reax_tools.reaction_snapshot.v1\",\n");
        fmt::print(fp, "  \"event_id\": {},\n", event.event_id);
        fmt::print(fp, "  \"frame_id\": {},\n", event.frame_id);
        fmt::print(fp, "  \"reactant_frame\": {},\n", event.frame_id - 1);
        fmt::print(fp, "  \"product_frame\": {},\n", event.frame_id);
        fmt::print(fp, "  \"time_fs\": {:.6f},\n", event.time_fs);
        fmt::print(fp, "  \"reaction_signature\": \"{}\",\n", json_escape(signature));
        fmt::print(fp, "  \"reaction_arity\": \"{}->{}\",\n", event.reactants.size(), event.products.size());
        fmt::print(fp, "  \"reactant_formulas\": \"{}\",\n", json_escape(event.get_reactants_string()));
        fmt::print(fp, "  \"product_formulas\": \"{}\",\n", json_escape(event.get_products_string()));
        fmt::print(fp, "  \"reactants\": [\n");
        print_molecules_json(fp, event.reactants, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"products\": [\n");
        print_molecules_json(fp, event.products, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"reactant_coords\": [\n");
        print_atoms_json(fp, event.reactants, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"product_coords\": [\n");
        print_atoms_json(fp, event.products, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"reactant_bonds\": [\n");
        print_bonds_json(fp, event.reactants, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"product_bonds\": [\n");
        print_bonds_json(fp, event.products, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"broken_bond_ids\": [\n");
        print_bond_change_ids_json(fp, reactant_bonds, product_bonds, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"created_bond_ids\": [\n");
        print_bond_change_ids_json(fp, product_bonds, reactant_bonds, "    ");
        fmt::print(fp, "\n  ],\n");
        fmt::print(fp, "  \"audit\": {{\n");
        fmt::print(fp, "    \"atom_ids_conserved\": {},\n", atom_ids_conserved ? "true" : "false");
        fmt::print(fp, "    \"reactant_atom_count\": {},\n", reactant_atoms.size());
        fmt::print(fp, "    \"product_atom_count\": {},\n", product_atoms.size());
        fmt::print(fp, "    \"reactant_bond_count\": {},\n", reactant_bonds.size());
        fmt::print(fp, "    \"product_bond_count\": {},\n", product_bonds.size());
        fmt::print(fp, "    \"broken_bond_count\": {},\n", broken_count);
        fmt::print(fp, "    \"created_bond_count\": {},\n", created_count);
        fmt::print(fp, "    \"allowed_arity\": true,\n");
        fmt::print(fp, "    \"signature_instance_index\": {}\n", signature_counts[signature]);
        fmt::print(fp, "  }}\n");
        fmt::print(fp, "}}\n");
        fclose(fp);

        fmt::print(manifest,
                   "{},{},\"{}\",{}->{},\"{}\",{},{},{},{},{}\n",
                   event.event_id,
                   event.frame_id,
                   signature,
                   event.reactants.size(),
                   event.products.size(),
                   rel_path,
                   atom_ids_conserved ? 1 : 0,
                   reactant_atoms.size(),
                   product_atoms.size(),
                   broken_count,
                   created_count);
        written++;
    }

    fclose(manifest);
    fmt::print("Saved {} reaction snapshot packages to {}\n", written, directory);
}

void ReactionTracker::save_raw_event_pairs(const std::string& filepath) {
    if (events.empty()) return;

    FILE* fp = create_file(filepath);

    fmt::print(fp, "event_id,frame,source_id,target_id,source_label,target_label,tracked_source_id,tracked_target_id,atom_overlap\n");

    int pair_count = 0;
    for (const auto& event : events) {
        for (const auto* reactant : event.reactants) {
            for (const auto* product : event.products) {
                int overlap = 0;
                const auto& smaller = reactant->atom_ids.size() < product->atom_ids.size()
                    ? reactant->atom_ids
                    : product->atom_ids;
                const auto& larger = reactant->atom_ids.size() < product->atom_ids.size()
                    ? product->atom_ids
                    : reactant->atom_ids;

                for (int atom_id : smaller) {
                    if (larger.count(atom_id)) {
                        overlap++;
                    }
                }

                if (overlap <= 0) continue;

                fmt::print(fp, "{},{},{},{},{},{},{},{},{}\n",
                           event.event_id,
                           event.frame_id,
                           molecule_identity_id(reactant->formula),
                           molecule_identity_id(product->formula),
                           reactant->formula,
                           product->formula,
                           reactant->id,
                           product->id,
                           overlap);
                pair_count++;
            }
        }
    }

    fclose(fp);
    fmt::print("Saved {} raw reaction event atom-overlap pairs to {}\n", pair_count, filepath);
}

void ReactionTracker::save_transfer_flow(const std::string& filepath) {
    if (events.empty()) return;

    struct TransferRecord {
        unsigned int source_id = 0;
        unsigned int target_id = 0;
        std::string source_label;
        std::string target_label;
        int count = 0;
        int atom_transfer = 0;
    };

    std::map<std::string, TransferRecord> records;

    for (const auto& event : events) {
        for (const auto* reactant : event.reactants) {
            for (const auto* product : event.products) {
                int overlap = 0;
                const auto& smaller = reactant->atom_ids.size() < product->atom_ids.size()
                    ? reactant->atom_ids
                    : product->atom_ids;
                const auto& larger = reactant->atom_ids.size() < product->atom_ids.size()
                    ? product->atom_ids
                    : reactant->atom_ids;

                for (int atom_id : smaller) {
                    if (larger.count(atom_id)) {
                        overlap++;
                    }
                }

                if (overlap <= 0) continue;

                unsigned int source_id = molecule_identity_id(reactant->formula);
                unsigned int target_id = molecule_identity_id(product->formula);
                std::string key = std::to_string(source_id) + " -> " + std::to_string(target_id);
                auto& record = records[key];
                record.source_id = source_id;
                record.target_id = target_id;
                record.source_label = reactant->formula;
                record.target_label = product->formula;
                record.count++;
                record.atom_transfer += overlap;
            }
        }
    }

    std::vector<TransferRecord> sorted_records;
    for (const auto& [key, record] : records) {
        sorted_records.push_back(record);
    }
    std::sort(sorted_records.begin(), sorted_records.end(),
        [](const auto& a, const auto& b) {
            if (a.source_id != b.source_id) return a.source_id < b.source_id;
            return a.target_id < b.target_id;
        });

    FILE* fp = create_file(filepath);
    fmt::print(fp, "source_id,target_id,source_label,target_label,count,atom_transfer,self_loop\n");
    for (const auto& record : sorted_records) {
        fmt::print(fp, "{},{},{},{},{},{},{}\n",
                   record.source_id,
                   record.target_id,
                   record.source_label,
                   record.target_label,
                   record.count,
                   record.atom_transfer,
                   record.source_id == record.target_id ? 1 : 0);
    }

    fclose(fp);
    fmt::print("Saved {} raw transfer-flow edges to {}\n", sorted_records.size(), filepath);
}

void ReactionTracker::save_events(const std::string& filepath) {
    if (events.empty()) return;
    
    // Count reaction frequencies and collect unique reactions
    struct ReactionRecord {
        std::string reactants_str;
        std::string products_str;
        int count;
        int first_frame;
        int last_frame;
    };
    
    std::map<std::string, ReactionRecord> reaction_map;
    
    for (const auto& event : events) {
        std::vector<std::string> r_species = split(event.get_reactants_string(), "+");
        std::vector<std::string> p_species = split(event.get_products_string(), "+");

        if (!cancel_common_species(r_species, p_species)) continue;

        if (!check_formula_atom_balance(r_species, p_species)) continue;

        std::sort(r_species.begin(), r_species.end());
        std::sort(p_species.begin(), p_species.end());

        std::string canon_reactants, canon_products;
        for (size_t i = 0; i < r_species.size(); ++i) {
            if (i > 0) canon_reactants += "+";
            canon_reactants += r_species[i];
        }
        for (size_t i = 0; i < p_species.size(); ++i) {
            if (i > 0) canon_products += "+";
            canon_products += p_species[i];
        }

        std::string key = canon_reactants + " -> " + canon_products;
        auto it = reaction_map.find(key);
        if (it == reaction_map.end()) {
            reaction_map[key] = {
                canon_reactants,
                canon_products,
                1,
                event.frame_id,
                event.frame_id
            };
        } else {
            it->second.count++;
            it->second.last_frame = event.frame_id;
        }
    }
    
    // Convert to vector and sort by frequency (descending)
    std::vector<ReactionRecord> sorted_reactions;
    for (const auto& [key, record] : reaction_map) {
        sorted_reactions.push_back(record);
    }
    std::sort(sorted_reactions.begin(), sorted_reactions.end(),
              [](const auto& a, const auto& b) { return a.count > b.count; });
    
    FILE* fp = create_file(filepath);
    
    // Header (without atom_transfer)
    fmt::print(fp, "rank,frequency,first_frame,last_frame,n_reactants,n_products,reactants,products\n");
    
    // Reactions sorted by frequency
    int rank = 1;
    for (const auto& rec : sorted_reactions) {
        // Count reactants and products
        int n_reactants = 1;
        int n_products = 1;
        for (char c : rec.reactants_str) if (c == '+') n_reactants++;
        for (char c : rec.products_str) if (c == '+') n_products++;
        
        fmt::print(fp, "{},{},{},{},{},{},\"{}\",\"{}\"\n",
                   rank++,
                   rec.count,
                   rec.first_frame,
                   rec.last_frame,
                   n_reactants,
                   n_products,
                   rec.reactants_str,
                   rec.products_str);
    }
    
    fclose(fp);
    fmt::print("Saved {} unique reaction types (sorted by frequency) to {}\n", sorted_reactions.size(), filepath);
}

void ReactionTracker::save_molecule_lifetimes(const std::string& filepath) {
    if (all_molecules.empty()) return;
    
    FILE* fp = create_file(filepath);
    
    // Header
    fmt::print(fp, "mol_id,formula,first_frame,last_frame,lifetime_frames,n_predecessors,n_successors\n");
    
    // Molecules
    for (const auto* mol : all_molecules) {
        int lifetime = mol->end_frame - mol->start_frame;
        fmt::print(fp, "{},{},{},{},{},{},{}\n",
                   mol->id,
                   mol->formula,
                   mol->start_frame,
                   mol->is_active() ? last_processed_frame : mol->end_frame,
                   lifetime,
                   mol->predecessors.size(),
                   mol->successors.size());
    }
    
    fclose(fp);
    fmt::print("Saved {} molecule lifetimes to {}\n", all_molecules.size(), filepath);
}

void ReactionTracker::brief_report() const {
    fmt::print("\n=== Reaction Tracker Report ===\n");
    fmt::print("Total molecules tracked: {}\n", all_molecules.size());
    fmt::print("Active molecules: {}\n", active_molecules.size());
    fmt::print("Reaction events detected: {}\n", events.size());
    
    if (!events.empty()) {
        fmt::print("\nTop 10 most common reactions:\n");
        
        // Count reaction types
        std::map<std::string, int> reaction_counts;
        for (const auto& event : events) {
            std::vector<std::string> r_species = split(event.get_reactants_string(), "+");
            std::vector<std::string> p_species = split(event.get_products_string(), "+");

            if (!cancel_common_species(r_species, p_species)) continue;
            if (!check_formula_atom_balance(r_species, p_species)) continue;

            std::sort(r_species.begin(), r_species.end());
            std::sort(p_species.begin(), p_species.end());

            std::string canon_reactants, canon_products;
            for (size_t i = 0; i < r_species.size(); ++i) {
                if (i > 0) canon_reactants += "+";
                canon_reactants += r_species[i];
            }
            for (size_t i = 0; i < p_species.size(); ++i) {
                if (i > 0) canon_products += "+";
                canon_products += p_species[i];
            }

            std::string key = canon_reactants + " -> " + canon_products;
            reaction_counts[key]++;
        }
        
        // Sort by count
        std::vector<std::pair<std::string, int>> sorted(reaction_counts.begin(), reaction_counts.end());
        std::sort(sorted.begin(), sorted.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        int max_display = std::min(10, (int)sorted.size());
        for (int i = 0; i < max_display; ++i) {
            fmt::print("  {}: {} times\n", sorted[i].first, sorted[i].second);
        }
    }
}
