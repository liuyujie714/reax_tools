#include "reax_flow.h"

#include <algorithm>
#include <fstream>
#include <mutex>
#include <unordered_map>

#include "argparser.h"
#include "fmt/format.h"
#include "string_tools.h"
#include "universe.h"

static std::string json_escape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (char c : value) {
        switch (c) {
        case '"': escaped += "\\\""; break;
        case '\\': escaped += "\\\\"; break;
        case '\n': escaped += "\\n"; break;
        case '\r': escaped += "\\r"; break;
        case '\t': escaped += "\\t"; break;
        default: escaped += c; break;
        }
    }
    return escaped;
}

/**
 * @brief Constructs a Node object representing a molecule in the reaction graph
 * @param mol Pointer to the molecule to be represented by this node
 * @throws ReaxFlowException if mol is nullptr
 * @note Creates a deep copy of the molecule to ensure the original can be
 * safely destroyed
 */
Node::Node(Molecule* mol) {
    // deep copy a molecule in system, thus the origin one can be detroyed safely.
    if (mol == nullptr) {
        throw ReaxFlowException("Invalid molecule when construct node.");
    }

    molecule = new Molecule(*mol);
    hash = molecule->hash;
}

/**
 * @brief Destructor for Node object
 * @note Safely deallocates the molecule pointer to prevent memory leaks
 */
Node::~Node() {
    if (molecule != nullptr) {
        delete molecule;
    }
}

/**
 * @brief Updates the degree statistics for this node based on reaction
 * participation
 * @param source_or_target True if this node is the source (outgoing), false if
 * target (incoming)
 * @param count Number of reactions involving this node
 * @param atom_transfer_count Total number of atoms transferred in reactions
 * involving this node
 * @note Updates topological degrees (edge counts) and weighted degrees (reaction counts)
 * @note degree = precursor_count + derivative_count (topological)
 * @note reaction_count = in_reaction_count + out_reaction_count (weighted)
 */
void Node::add_degrees(bool source_or_target, unsigned int count, unsigned int atom_transfer_count) {
    if (source_or_target) {
        // As reactant/source: forms derivatives
        derivative_count++;                     // topological
        derivative_reactions += count;          // weighted
        derivative_atom_transfer += atom_transfer_count;
    }
    else {
        // As product/target: has precursors
        precursor_count++;                      // topological
        precursor_reactions += count;           // weighted
        precursor_atom_transfer += atom_transfer_count;
    }

    // Update totals
    degree = precursor_count + derivative_count;
    reaction_count = precursor_reactions + derivative_reactions;
    atom_transfer = precursor_atom_transfer + derivative_atom_transfer;
}

/**
 * @brief Constructs an Edge object representing a reaction between two
 * molecules
 * @param from_node Pointer to the source node (reactant)
 * @param to_node Pointer to the target node (product)
 * @throws ReaxFlowException if either node is nullptr
 * @note Automatically generates a hash for efficient lookup
 */
Edge::Edge(Node* from_node, Node* to_node) {
    if (from_node == nullptr || to_node == nullptr) {
        throw ReaxFlowException("Invalid node when construct edge.");
    }

    source = from_node;
    target = to_node;
    hash = get_edge_hash(source, target);
}

/**
 * @brief Destructor for Edge object
 * @note No dynamic memory to deallocate, but provides virtual destructor for
 * inheritance
 */
Edge::~Edge() {}

/**
 * @brief Default constructor for ReaxFlow object
 * @note Initializes empty reaction graph with thread-safe mutex
 */
ReaxFlow::ReaxFlow() {}

/**
 * @brief Destructor for ReaxFlow object
 * @note Safely deallocates all nodes and edges to prevent memory leaks
 */
ReaxFlow::~ReaxFlow() {
    for (auto& edge : edges) {
        if (edge != nullptr) {
            delete edge;
        }
    }
    edges.clear();

    for (auto& node : nodes) {
        if (node != nullptr) {
            delete node;
        }
    }
    nodes.clear();

    // Clear the hash map
    molecule_hash_to_node.clear();
}

/**
 * @brief Retrieves a node by molecule pointer using hash-based lookup
 * @param mol Pointer to the molecule to search for
 * @return Node* Pointer to the node if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Node* ReaxFlow::get_node(Molecule* mol) {
    if (mol == nullptr) {
        return nullptr;
    }

    // Use hash map for O(1) lookup
    auto it = molecule_hash_to_node.find(mol->hash);
    if (it != molecule_hash_to_node.end()) {
        return it->second;
    }

    return nullptr;
}

/**
 * @brief Retrieves a node by hash value using hash-based lookup
 * @param hash Hash value of the molecule to search for
 * @return Node* Pointer to the node if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Node* ReaxFlow::get_node(unsigned int hash) {
    auto it = molecule_hash_to_node.find(hash);
    if (it != molecule_hash_to_node.end()) {
        return it->second;
    }
    return nullptr;
}

/**
 * @brief Retrieves an edge between two nodes using hash-based lookup
 * @param from_node Source node of the edge
 * @param to_node Target node of the edge
 * @return Edge* Pointer to the edge if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Edge* ReaxFlow::get_edge(Node* from_node, Node* to_node) {
    if (from_node == nullptr || to_node == nullptr) {
        return nullptr;
    }

    unsigned int hash = get_edge_hash(from_node, to_node);
    auto it = edge_hash_to_edge.find(hash);
    if (it != edge_hash_to_edge.end()) {
        return it->second;
    }
    return nullptr;
}

/**
 * @brief Retrieves an edge by hash value using hash-based lookup
 * @param hash Hash value of the edge to search for
 * @return Edge* Pointer to the edge if found, nullptr otherwise
 * @note Uses O(1) hash map lookup for efficient retrieval
 */
Edge* ReaxFlow::get_edge(unsigned int hash) {
    auto it = edge_hash_to_edge.find(hash);
    if (it != edge_hash_to_edge.end()) {
        return it->second;
    }
    return nullptr;
}

/**
 * @brief Adds a molecule to the reaction flow system
 * @param mol Pointer to the molecule to be added
 * @return Node* Pointer to the node representing the molecule
 * @note If the molecule already exists, returns its existing node. Otherwise
 * creates a new node.
 * @note Thread-safe implementation with automatic node creation
 */
Node* ReaxFlow::add_molecule(Molecule* mol) {
    if (mol == nullptr) {
        return nullptr;
    }

    // Check if molecule already exists using hash map
    auto it = molecule_hash_to_node.find(mol->hash);
    if (it != molecule_hash_to_node.end()) {
        return it->second;  // Return existing node
    }

    // Create new node
    Node* new_node = new Node(mol);
    nodes.insert(new_node);
    molecule_hash_to_node[mol->hash] = new_node;

    return new_node;
}

Node* ReaxFlow::register_molecule(Molecule* mol) {
    std::lock_guard<std::mutex> lock(reaxflow_mutex);
    return add_molecule(mol);
}

/**
 * @brief Adds a reaction between two molecules to the system
 * @param frame_id Frame number of the reaction (currently unused but preserved
 * for future use)
 * @param atom_transfer_count Number of atoms transferred in the reaction
 * @param source Source molecule of the reaction (reactant)
 * @param target Target molecule of the reaction (product)
 * @return Edge* Pointer to the edge representing the reaction
 * @note Thread-safe implementation with mutex lock
 * @note If the reaction already exists, updates the existing edge counts
 * @note Automatically creates nodes for molecules if they don't exist
 */
Edge* ReaxFlow::add_reaction(const int& frame_id, const int& atom_transfer_count, Molecule* source, Molecule* target) {
    // lock for parallel, dynamic find & create in unordered_set / unordered_map
    // is thread unsafe.
    std::lock_guard<std::mutex> lock(reaxflow_mutex);

    Node* source_node = add_molecule(source);
    Node* target_node = add_molecule(target);

    if (source_node == nullptr || target_node == nullptr) {
        return nullptr;
    }

    // Check if edge already exists
    Edge* existing_edge = get_edge(source_node, target_node);
    if (existing_edge != nullptr) {
        // Update existing edge
        existing_edge->count++;
        existing_edge->atom_transfer += atom_transfer_count;

        return existing_edge;
    }

    // Create new edge
    Edge* new_edge = new Edge(source_node, target_node);
    new_edge->count++;
    new_edge->atom_transfer += atom_transfer_count;

    edges.insert(new_edge);
    edge_hash_to_edge[new_edge->hash] = new_edge;

    return new_edge;
}

/**
 * @brief Sorts nodes and edges by their degree/count for efficient reporting
 * @note Calls calc_node_degrees() to ensure degree data is current
 * @note Clears and rebuilds sorted_nodes and sorted_edges vectors
 * @note Sorts in descending order (highest degree/count first)
 */
void ReaxFlow::update_graph() {
    // Clear all degree statistics and adjacency relationships
    for (auto& node : nodes) {
        // Topological degrees (precursors and derivatives)
        node->degree = 0;
        node->precursor_count = 0;
        node->derivative_count = 0;

        // Weighted degrees (reaction counts)
        node->reaction_count = 0;
        node->precursor_reactions = 0;
        node->derivative_reactions = 0;

        // Atom transfers
        node->atom_transfer = 0;
        node->precursor_atom_transfer = 0;
        node->derivative_atom_transfer = 0;

        // Adjacency (precursors and derivatives)
        node->from_nodes.clear();
        node->to_nodes.clear();
    }

    // Build adjacency relationships and calculate degrees
    for (const auto& edge : edges) {
        edge->source->add_degrees(true, edge->count, edge->atom_transfer);
        edge->target->add_degrees(false, edge->count, edge->atom_transfer);

        // Build adjacency relationships
        edge->source->to_nodes.insert(edge->target);
        edge->target->from_nodes.insert(edge->source);
    }

    sorted_nodes.clear();
    sorted_edges.clear();

    for (const auto& node : nodes) {
        sorted_nodes.emplace_back(std::pair(node, node->degree));
    }
    std::sort(sorted_nodes.begin(), sorted_nodes.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });

    for (const auto& edge : edges) {
        sorted_edges.emplace_back(std::pair(edge, edge->count));
    }
    std::sort(sorted_edges.begin(), sorted_edges.end(),
        [](const auto& a, const auto& b) { return a.second > b.second; });
}

/**
 * @brief Generates a brief report of the reaction flow statistics
 * @note Displays top 10 molecules by degree and top 20 reactions by count
 * @note Shows in-degree and out-degree for molecules, reaction counts and atom
 * transfers for edges
 * @note Automatically calls sort_nodes_and_edges() to ensure data is current
 */
void ReaxFlow::brief_report() {
    update_graph();

    unsigned int max_node_display = std::min(10, int(sorted_nodes.size()));
    unsigned int max_edge_display = std::min(20, int(sorted_edges.size()));

    fmt::print("\n=== Reaction Flow Report ===\n");
    fmt::print("Top {} key molecules:\n", max_node_display);
    fmt::print("{:<12s}{:<12s}{:<12s}\n", "molecule", "precursors", "derivatives");
    Node* tmp_node = nullptr;
    for (size_t i = 0; i < max_node_display; i++) {
        tmp_node = sorted_nodes[i].first;
        fmt::print("{:<12s}{:<12d}{:<12d}\n", tmp_node->molecule->formula, tmp_node->precursor_count, tmp_node->derivative_count);
    }

    fmt::print("\n");
    fmt::print("Top {} key flows:\n", max_edge_display);
    fmt::print("{:<12s}{:<12s}{:<8s}{:<15s}\n", "from", "to", "count", "atom transfered");
    Node* tmp_source = nullptr;
    Node* tmp_target = nullptr;
    for (size_t i = 0; i < max_edge_display; i++) {
        tmp_source = sorted_edges[i].first->source;
        tmp_target = sorted_edges[i].first->target;
        fmt::print("{} -> {} = R:{} AT:{}\n", tmp_source->molecule->formula, tmp_target->molecule->formula,
            sorted_edges[i].second, sorted_edges[i].first->atom_transfer);
    }

    tmp_node = nullptr;
    tmp_source = nullptr;
    tmp_target = nullptr;
}

/**
 * @brief Writes the reaction flow graph to a DOT file for visualization
 * @param basename Path to the output DOT file
 * @param edges_to_write Set of edges to include in the graph
 * @param write_atom_transfer Whether to use atom transfer counts as edge
 * weights
 * @param layout Graph layout algorithm to use (e.g., "circo", "dot", "neato")
 * @note Creates a directed graph with nodes representing molecules and edges
 * representing reactions
 * @note Edge thickness is proportional to reaction count (logarithmic scale)
 * @note Highlights top 25% of edges in goldenrod color
 */
void ReaxFlow::save_molecules_json() {
    update_graph();

    std::vector<Node*> sorted_nodes_vec(nodes.begin(), nodes.end());
    std::sort(sorted_nodes_vec.begin(), sorted_nodes_vec.end(),
        [](Node* a, Node* b) {
            return a->hash < b->hash;
        });

    FILE* fp = create_file("molecules.json");
    fmt::print(fp, "{{\n");
    fmt::print(fp, "  \"schema_version\": \"1.0\",\n");
    fmt::print(fp, "  \"identity_model\": \"formula-hash-v1\",\n");
    fmt::print(fp, "  \"molecules\": [\n");

    for (size_t node_i = 0; node_i < sorted_nodes_vec.size(); ++node_i) {
        Node* node = sorted_nodes_vec[node_i];
        Molecule* mol = node->molecule;

        std::vector<Atom*> atoms(mol->mol_atoms.begin(), mol->mol_atoms.end());
        std::sort(atoms.begin(), atoms.end(), [](Atom* a, Atom* b) { return a->id < b->id; });

        std::vector<Bond*> bonds(mol->mol_bonds.begin(), mol->mol_bonds.end());
        std::sort(bonds.begin(), bonds.end(), [](Bond* a, Bond* b) {
            int a_min = std::min(a->atom_i->id, a->atom_j->id);
            int a_max = std::max(a->atom_i->id, a->atom_j->id);
            int b_min = std::min(b->atom_i->id, b->atom_j->id);
            int b_max = std::max(b->atom_i->id, b->atom_j->id);
            if (a_min != b_min) return a_min < b_min;
            return a_max < b_max;
        });

        fmt::print(fp, "    {{\n");
        fmt::print(fp, "      \"id\": \"{}\",\n", node->hash);
        fmt::print(fp, "      \"formula\": \"{}\",\n", json_escape(mol->formula));
        fmt::print(fp, "      \"atom_counts\": {{");
        bool first_count = true;
        for (const auto& [element, count] : mol->types_nums) {
            if (!first_count) fmt::print(fp, ", ");
            fmt::print(fp, "\"{}\": {}", json_escape(element), count);
            first_count = false;
        }
        fmt::print(fp, "}},\n");

        fmt::print(fp, "      \"stats\": {{\n");
        fmt::print(fp, "        \"total_connections\": {},\n", node->degree);
        fmt::print(fp, "        \"precursors\": {},\n", node->precursor_count);
        fmt::print(fp, "        \"derivatives\": {}\n", node->derivative_count);
        fmt::print(fp, "      }},\n");

        fmt::print(fp, "      \"example\": {{\n");
        fmt::print(fp, "        \"atoms\": [\n");
        for (size_t i = 0; i < atoms.size(); ++i) {
            Atom* atom = atoms[i];
            fmt::print(fp, "          {{\"id\": {}, \"element\": \"{}\"}}{}{}\n",
                atom->id,
                json_escape(atom->type_name),
                i + 1 < atoms.size() ? "," : "",
                "");
        }
        fmt::print(fp, "        ],\n");

        fmt::print(fp, "        \"bonds\": [\n");
        for (size_t i = 0; i < bonds.size(); ++i) {
            Bond* bond = bonds[i];
            fmt::print(fp, "          {{\"a\": {}, \"b\": {}, \"order\": {}}}{}\n",
                bond->atom_i->id,
                bond->atom_j->id,
                bond->order,
                i + 1 < bonds.size() ? "," : "");
        }
        fmt::print(fp, "        ]\n");
        fmt::print(fp, "      }}\n");
        fmt::print(fp, "    }}{}\n", node_i + 1 < sorted_nodes_vec.size() ? "," : "");
    }

    fmt::print(fp, "  ]\n");
    fmt::print(fp, "}}\n");
    fclose(fp);
}

/**
 * @brief Merges molecules by element count ranges for simplified analysis
 * @param target_element Element symbol to group by (e.g., "C", "H", "O")
 * @param ranges Vector of range boundaries for grouping (e.g., {0, 5, 10}
 * creates groups 0-4, 5-9, 10+)
 * @note Creates grouped molecules with formula prefix "grp_" followed by
 * element and range
 * @note Molecules with target element count in each range are merged into a
 * single node
 * @note Useful for simplifying complex reaction networks by grouping similar
 * molecules
 */
void ReaxFlow::merge_by_element(std::string target_element, std::vector<int> ranges) {
    std::unordered_set<std::string> all_formulas;
    std::unordered_map<std::string, std::unordered_set<std::string>> formulas_map;
    std::map<std::string, int> elements_weights;

    for (const auto& node : nodes) {
        all_formulas.insert(node->molecule->formula);
    }

    for (size_t i = 0; i < ranges.size(); i++) {
        int start;
        int end;
        std::string new_formula;

        if (i < ranges.size() - 1) {
            start = ranges[i];
            end = ranges[i + 1] - 1;
            new_formula = fmt::format("grp_{}{}-{}", target_element, start, end);
        }
        else {
            start = ranges[i];
            end = 100000;  // If a molecule have > 10000 atoms, that's user's bad input,
            // will not provide anything makes sense from the beginning.
            new_formula = fmt::format("grp_{}{}-max", target_element, start);
        }

        formulas_map[new_formula] = {};

        for (const auto& formula : all_formulas) {
            if (starts_with("grp_", formula)) {
                continue;
            }

            elements_weights = parse_formula(formula);
            for (const auto& [elem, weight] : elements_weights) {
                if (target_element == elem && weight >= start && weight <= end) {
                    formulas_map[new_formula].insert(formula);
                }
            }
        }
    }

    for (const auto& [new_formula, old_formulas_set] : formulas_map) {
        merge_formulas(old_formulas_set, new_formula);
    }
}

/**
 * @brief Merges multiple molecules into a single representative node
 * @param formulas_set Set of formula strings to merge
 * @param new_formula New formula string for the merged node
 * @note Combines all nodes with formulas in formulas_set into a single node
 * @note Updates all edges to connect to the new merged node
 * @note Preserves reaction counts and atom transfer data
 * @note Deallocates old nodes and edges to prevent memory leaks
 */
void ReaxFlow::merge_formulas(const std::unordered_set<std::string>& formulas_set, const std::string& new_formula) {
    // Get target nodes to merge
    std::vector<Node*> nodes_to_merge;
    for (const auto& node : nodes) {
        if (formulas_set.count(node->molecule->formula)) {
            nodes_to_merge.push_back(node);
        }
    }

    // Nothing to merge
    if (nodes_to_merge.size() < 1) return;

    // Create a new merged node
    Node* merged_node = new Node(nodes_to_merge[0]->molecule);
    merged_node->molecule->formula = new_formula;
    merged_node->hash = std::hash<std::string>()(new_formula);

    // Add the merged node
    nodes.insert(merged_node);
    molecule_hash_to_node[merged_node->hash] = merged_node;

    // Collect all edges involving the nodes to be merged
    std::vector<Edge*> edges_to_remove;
    std::vector<Edge*> edges_to_add;

    for (const auto& node_to_merge : nodes_to_merge) {
        // Find all edges involving this node
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            Edge* edge = *it;

            if (edge->source == node_to_merge || edge->target == node_to_merge) {
                edges_to_remove.push_back(edge);

                // Create new edge with merged node
                Node* other_node = (edge->source == node_to_merge) ? edge->target : edge->source;
                Edge* new_edge = new Edge((edge->source == node_to_merge) ? merged_node : other_node,
                    (edge->target == node_to_merge) ? merged_node : other_node);
                new_edge->count = edge->count;
                new_edge->atom_transfer = edge->atom_transfer;
                edges_to_add.push_back(new_edge);
            }
        }

        // Remove the node
        nodes.erase(node_to_merge);
        molecule_hash_to_node.erase(node_to_merge->hash);
        delete node_to_merge;
    }

    // Remove old edges
    for (Edge* edge : edges_to_remove) {
        edges.erase(edges.find(edge));
        delete edge;
    }

    // Add new edges
    for (Edge* edge : edges_to_add) {
        edges.insert(edge);
    }
}
