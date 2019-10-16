#include "LatticeMesh2D.h"
#include "LatticeMesh3D.h"

int main(int argc, char *argv[])
{
    #if 0
    LatticeMesh<float, 2> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.05;
    simulationMesh.m_nFrames = 30;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();
    simulationMesh.initializeDeformation();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();
#else
     LatticeMesh<float, 3> simulationMesh;
    simulationMesh.m_cellSize = { 10, 10, 10 };
    simulationMesh.m_gridDX = 0.05;
    simulationMesh.m_nFrames = 30;
    simulationMesh.m_subSteps = 1;
    simulationMesh.m_frameDt = 0.1;

    // Initialize the simulation example
    simulationMesh.initialize();
    simulationMesh.initializeDeformation();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();
#endif
    return 0;
}

