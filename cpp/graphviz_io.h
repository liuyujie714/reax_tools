#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "fmt/core.h"
#include <graphviz/gvc.h>
#include <graphviz/cgraph.h>
#include <graphviz/gvplugin.h>

#ifdef __cplusplus
extern "C" {
#endif

extern gvplugin_library_t gvplugin_dot_layout_LTX_library;
extern gvplugin_library_t gvplugin_gdiplus_LTX_library;


/* \brief Convert Dot file(Graphviz) to PNG file with dpi
* \param dot file name
* \param png file name
* \param dpi DPI
*/
void dot_to_png(const std::string& dotfname, const std::string& pngfname, int dpi = 300)
{
    auto dot_file = (std::filesystem::path(OUTPUT_DIR) / dotfname).string();
    auto png_file = (std::filesystem::path(OUTPUT_DIR) / pngfname).string();

    // 1. Read the DOT file
	std::ifstream dot(dot_file);
	if (!dot.is_open())
	{
        fmt::print("Failed to open file: {}", dot_file);
        return;
	}
	std::stringstream buffer;
	buffer << dot.rdbuf();
	std::string dot_str = buffer.str();

    // 2. Init gvc lib & Parse the DOT string
    GVC_t* gvc = gvContext();
    gvAddLibrary(gvc, &gvplugin_dot_layout_LTX_library);
    gvAddLibrary(gvc, &gvplugin_gdiplus_LTX_library);
    Agraph_t* g = agmemread(const_cast<char*>(dot_str.c_str()));

    // 3. Set dpi & Layout the graph
    agattr(g, AGNODE, (char*)"fontname", "Arial");
    agattr(g, AGNODE, (char*)"fontsize", "14");
    agattr(g, AGRAPH, (char*)"dpi", std::to_string(dpi).c_str());
    gvLayout(gvc, g, "dot");

    // 4. Render to PNG
    gvRenderFilename(gvc, g, "png", png_file.c_str());

    // 5. Clean up memory
    gvFreeLayout(gvc, g);
    agclose(g);
    gvFreeContext(gvc);

    fmt::print("\n{} rendered to: {}\n", dot_file, png_file);
}

#ifdef __cplusplus
}
#endif
