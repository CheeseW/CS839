#pragma once

#include "FiniteElementMesh3D.h"

#include <Eigen/Dense>

template<class T, int d>
struct LatticeMesh;

template<class T>
struct LatticeMesh<T, 3> : public FiniteElementMesh<T,3>
{
    static constexpr int d = 3;
    using Base = FiniteElementMesh<T,d>;

    // from AnimatedMesh
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using VectorType = typename Base::VectorType;

    // from FiniteElementMesh
    using Base::initializeUndeformedConfiguration;
    using Base::m_stepEndTime;

    std::array<int, d> m_cellSize; // dimensions in grid cells
    T m_gridDX;

    const int m_pinchRadius;

    std::vector<VectorType> m_particleUndeformedX;
    std::vector<int> m_leftHandleIndices;
    std::vector<int> m_rightHandleIndices;
    VectorType m_leftHandleDisplacement;
    VectorType m_rightHandleDisplacement;
    
    LatticeMesh()
        :Base(1.e2, 1., 4., .05), m_pinchRadius(1)
    { 
        // TODO: check this 3D adaption is correct
        m_leftHandleDisplacement  = VectorType(-.2, 0., 0.);
        m_rightHandleDisplacement = VectorType( .2, 0., 0.);
    }

    void initialize()
    {
        initializeUSD("HW02.usda");
#if 0
        // Create a Cartesian lattice topology
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++){
            m_meshElements.emplace_back(
                std::array<int, 3>{
                    gridToParticleID(cell_i  , cell_j  ), 
                    gridToParticleID(cell_i+1, cell_j  ),
                    gridToParticleID(cell_i+1, cell_j+1)
                }
            );
            m_meshElements.emplace_back(
                std::array<int, 3>{
                    gridToParticleID(cell_i  , cell_j  ), 
                    gridToParticleID(cell_i+1, cell_j+1),
                    gridToParticleID(cell_i  , cell_j+1)
                }
            );
        }
#endif

     std::vector<std::array<int, 3>> activeCells;
    //std::map<std::array<int, 3>, int> activeNodes;

        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
        for(int cell_k = 0; cell_k < m_cellSize[1]; cell_k++){
                activeCells.push_back(std::array<int, 3>{cell_i, cell_j, cell_k});
        }

        for(const auto& cell: activeCells){
            std::array<int, 3> node;
            for(node[0] = cell[0]; node[0] <= cell[0]+1; node[0]++)
            for(node[1] = cell[1]; node[1] <= cell[1]+1; node[1]++)
            for(node[2] = cell[2]; node[2] <= cell[2]+1; node[2]++){
                auto search = m_activeNodes.find(node);
                if(search == m_activeNodes.end()){ // Particle not yet created at this lattice node location -> make one
                    m_activeNodes.insert({node, m_particleX.size()});
                    m_particleX.emplace_back(m_gridDX * T(node[0]), m_gridDX * T(node[1]), m_gridDX * T(node[2]));
                }
            }
        }

        for(const auto& cell: activeCells){
            int vertexIndices[2][2][2];
            for(int i = 0; i <= 1; i++)
            for(int j = 0; j <= 1; j++)
            for(int k = 0; k <= 1; k++){
                std::array<int, 3> node{cell[0] + i, cell[1] + j, cell[2] + k};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    vertexIndices[i][j][k] = search->second;
                else
                    throw std::logic_error("particle at cell vertex not found");
            }

            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][0], vertexIndices[1][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][0], vertexIndices[1][1][1], vertexIndices[1][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][0][1], vertexIndices[1][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][1], vertexIndices[0][0][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][1], vertexIndices[0][1][0], vertexIndices[0][1][1]});
            m_meshElements.push_back(std::array<int, 4>{ vertexIndices[0][0][0], vertexIndices[1][1][0], vertexIndices[0][1][0], vertexIndices[1][1][1]});
        }

        initializeTopology();
#if 0
        // Also initialize the associated particles
        for(int node_i = 0; node_i <= m_cellSize[0]; node_i++)
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++)
            m_particleX.emplace_back(m_gridDX * (T)node_i, m_gridDX * (T)node_j);
#endif
        initializeParticles();

        // Check particle indexing in mesh
        for(const auto& element: m_meshElements)
            for(const auto vertex: element)
                if(vertex < 0 || vertex >= m_particleX.size())
                    throw std::logic_error("mismatch between mesh vertex and particle array");

        // Initialize rest shape matrices and particle mass
        initializeUndeformedConfiguration();

        // Also record rest shape
        m_particleUndeformedX = m_particleX;

        // Identify particles on left and right handles
        // Assumed that the topology is a rectangular grid 
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++)
        for(int node_k = 0; node_k <= m_cellSize[2]; node_k++){
            {
             std::array<int, 3> node{0, node_j, node_k};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    m_leftHandleIndices.push_back(search->second);
                else
                    throw std::logic_error("particle at cell vertex not found");
                    
            }
             {
             std::array<int, 3> node{m_cellSize[0], node_j, node_k};
                auto search = m_activeNodes.find(node);
                if(search != m_activeNodes.end())
                    m_rightHandleIndices.push_back(search->second);
                else
                    throw std::logic_error("particle at cell vertex not found");
                    
            }

        }

    }

    void initializeDeformation()
    {
        // No need to apply any deformation; this example is driven by moving handles
    }

    void clearConstrainedParticles(std::vector<VectorType>& x) override
    { 
        for(const auto v: m_leftHandleIndices)
            x[v] = VectorType::Zero();
        for(const auto v: m_rightHandleIndices)
            x[v] = VectorType::Zero();       
    }

    void setBoundaryConditions() override
    { 
        T effectiveTime = std::min<T>(m_stepEndTime, 1.0);
        
        for(const auto v: m_leftHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_leftHandleDisplacement;
        }
        for(const auto v: m_rightHandleIndices){
            m_particleX[v] = m_particleUndeformedX[v] + effectiveTime * m_rightHandleDisplacement;
        }
    }

private:
   // inline int gridToParticleID(const int i, const int j) const { return i * (m_cellSize[1]+1) + j; }
   std::map<std::array<int, 3>, int> m_activeNodes;
};