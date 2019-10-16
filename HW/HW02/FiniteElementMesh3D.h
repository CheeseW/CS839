#pragma once

#include "AnimatedTetrahedronMesh.h"

#include "CGVectorWrapper.h"
#include "CGSystemWrapper.h"
#include "ConjugateGradient.h"

#include <Eigen/Dense>

//#define USE_LINEAR_ELASTICITY
#define USE_ST_VENANT_KIRCHHOFF
//#define USE_COROTATED_ELASTICITY
template<class T, int d>
struct FiniteElementMesh;

template<class T>
struct FiniteElementMesh<T,3> : public AnimatedTetrahedronMesh<T>
{
    using Base = AnimatedTetrahedronMesh<T>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using VectorType = typename Base::Vector3;
     static constexpr int d = 3;
    using MatrixType = Eigen::Matrix< T , d , d>;
   

    int m_nFrames;
    int m_subSteps;
    T m_frameDt;
    T m_stepDt;
    T m_stepEndTime;

    const T m_density;
    const T m_mu;
    const T m_lambda;
    const T m_rayleighCoefficient;
    const T m_singularValueThreshold;
    
    std::vector<T> m_particleMass;
    std::vector<MatrixType> m_DmInverse;
    std::vector<T> m_restVolume;
    
    FiniteElementMesh(const T density, const T mu, const T lambda, const T rayleighCoefficient)
        :m_density(density), m_mu(mu), m_lambda(lambda), m_rayleighCoefficient(rayleighCoefficient), m_singularValueThreshold(.2f)
    {}

    void initializeUndeformedConfiguration()
    {
        // Initialize rest shape and particle mass (based on constant density)
        m_particleMass.resize(m_particleX.size(), T()); // Initialize all particle masses to zero
        for(const auto& element: m_meshElements)
        {
            MatrixType Dm;
            for(int j = 0; j < d; j++)
                Dm.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            T restVolume = .5 * Dm.determinant();
            if(restVolume < 0)
                throw std::logic_error("Inverted element");
            m_DmInverse.emplace_back(Dm.inverse());
            m_restVolume.push_back(restVolume);
            T elementMass = m_density * restVolume;
            for(const int v: element)
                m_particleMass[v] += (1./3.) * elementMass; // TODO: make scale relates to d
        }
    }
    
    void addElasticForce(std::vector<VectorType>& f) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Compute deformation gradient
            MatrixType Ds;
            for(int j = 0; j < d; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            MatrixType F = Ds * m_DmInverse[e];

            // Compute SVD
            Eigen::JacobiSVD<MatrixType> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            MatrixType U = svd.matrixU();
            MatrixType V = svd.matrixV();
            VectorType vSigma = svd.singularValues();
            if ( U.determinant() < 0. ) {
                if ( V.determinant() < 0. ) {
                    // Both determinants negative, just negate 2nd column on both
                    U.col(d-1) *= -1.f;
                    V.col(d-1) *= -1.f;
                }
                else {
                    // Only U has negative determinant, negate 2nd column and second singular value
                    U.col(d-1) *= -1.f;
                    vSigma[d-1] = -vSigma[d-1];
                }
            }
            else
                if ( V.determinant() < 0.) {
                    // Only V has negative determinant, negate 2nd column and second singular value
                    V.col(d-1) *= -1.f;
                    vSigma[d-1] = -vSigma[d-1];
                }
            if ( (F-U*vSigma.asDiagonal()*V.transpose()).norm() > 1e-5 )
                throw std::logic_error("SVD error");

            // Apply thresholding of singular values, and re-constitute F
            for (int v = 0; v < d; v++)
                vSigma[v] = std::max<T>(m_singularValueThreshold, vSigma[v]);
            MatrixType Sigma = vSigma.asDiagonal();
            F = U * Sigma * V.transpose();
            
#ifdef USE_LINEAR_ELASTICITY
            MatrixType strain = .5 * (F + F.transpose()) - MatrixType::Identity();
            MatrixType P = 2. * m_mu * strain + m_lambda * strain.trace() * MatrixType::Identity();
#endif

#ifdef USE_ST_VENANT_KIRCHHOFF
            MatrixType E = .5 * ( F.transpose() * F - MatrixType::Identity());
            MatrixType P = F * (2. * m_mu * E + m_lambda * E.trace() * MatrixType::Identity());
#endif

#ifdef USE_COROTATED_ELASTICITY
            VectorType vStrain = vSigma - VectorType::Ones();
            VectorType vP = 2. * m_mu * vStrain + m_lambda * vStrain.sum() * VectorType::Ones();
            MatrixType P = U * vP.asDiagonal() * V.transpose();
#endif

            MatrixType H = -m_restVolume[e] * P * m_DmInverse[e].transpose();
            
            for(int j = 0; j < d; j++){
                f[element[j+1]] += H.col(j);
                f[element[0]] -= H.col(j);
            }
        }
    }

    void addProductWithStiffnessMatrix(std::vector<VectorType>& dx, std::vector<VectorType>& df, const T scale) const
    {
        for(int e = 0; e < m_meshElements.size(); e++)
        {
            const auto& element = m_meshElements[e];

            // Compute deformation gradient
            MatrixType Ds;
            for(int j = 0; j < d; j++)
                Ds.col(j) = m_particleX[element[j+1]]-m_particleX[element[0]];
            MatrixType F = Ds * m_DmInverse[e];

            // Compute SVD
            Eigen::JacobiSVD<MatrixType> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
            MatrixType U = svd.matrixU();
            MatrixType V = svd.matrixV();
            VectorType vSigma = svd.singularValues();
            if ( U.determinant() < 0. ) {
                if ( V.determinant() < 0. ) {
                    // Both determinants negative, just negate 2nd column on both
                    U.col(d-1) *= -1.f;
                    V.col(d-1) *= -1.f;
                }
                else {
                    // Only U has negative determinant, negate 2nd column and second singular value
                    U.col(d-1) *= -1.f;
                    vSigma[d-1] = -vSigma[d-1];
                }
            }
            else
                if ( V.determinant() < 0.) {
                    // Only V has negative determinant, negate 2nd column and second singular value
                    V.col(d-1) *= -1.f;
                    vSigma[d-1] = -vSigma[d-1];
                }
            if ( (F-U*vSigma.asDiagonal()*V.transpose()).norm() > 1e-5 )
                throw std::logic_error("SVD error");

            // Apply thresholding of singular values, and re-constitute F
            for (int v = 0; v < d; v++)
                vSigma[v] = std::max<T>(m_singularValueThreshold, vSigma[v]);
            MatrixType Sigma = vSigma.asDiagonal();
            F = U * Sigma * V.transpose();
            
            // Compute differential(s)
            MatrixType dDs;
            for(int j = 0; j < d; j++)
                dDs.col(j) = dx[element[j+1]]-dx[element[0]];
            MatrixType dF = dDs * m_DmInverse[e];

#ifdef USE_LINEAR_ELASTICITY
            MatrixType dstrain = .5 * (dF + dF.transpose());
            MatrixType dP = scale * (2. * m_mu * dstrain + m_lambda * dstrain.trace() * MatrixType::Identity());
#endif

#ifdef USE_ST_VENANT_KIRCHHOFF
            MatrixType E = .5 * ( F.transpose() * F - MatrixType::Identity());
            MatrixType dE = .5 * ( dF.transpose() * F + F.transpose() * dF);
            MatrixType dP = dF * (2. * m_mu *  E + m_lambda *  E.trace() * MatrixType::Identity()) +
                           F * (2. * m_mu * dE + m_lambda * dE.trace() * MatrixType::Identity());

#endif

#ifdef USE_COROTATED_ELASTICITY // TODO: make this part work for 3d
            // Construct diagonalized dP/dF tensor
            MatrixType A, B12;
            A(0, 0) = A(1, 1) = 2. * m_mu + m_lambda;
            A(0, 1) = A(1, 0) = m_lambda;
            VectorType vStrain = vSigma - VectorType::Ones();
            T q = std::max<T>( m_lambda * vStrain.sum() - 2. * m_mu, 0. ); // Positive definiteness fix
            B12(0, 0) = B12(1, 1) = m_mu + q;
            B12(0, 1) = B12(1, 0) = m_mu - q;

            // Apply tensor (with required rotations)
            MatrixType dF_hat = U.transpose() * dF * V;

            VectorType vdP_A = A * VectorType( dF_hat(0, 0), dF_hat(1, 1) );
            VectorType vdP_B12 = B12 * VectorType( dF_hat(0, 1), dF_hat(1, 0) );
            MatrixType dP_hat;
            dP_hat(0, 0) = vdP_A[0];
            dP_hat(1, 1) = vdP_A[1];
            dP_hat(0, 1) = vdP_B12[0];
            dP_hat(1, 0) = vdP_B12[1];

            MatrixType dP = U * dP_hat * V.transpose();
#endif
            
            MatrixType dH = m_restVolume[e] * dP * m_DmInverse[e].transpose();
            
            for(int j = 0; j < d; j++){
                df[element[j+1]] += dH.col(j);
                df[element[0]] -= dH.col(j);
            }
        }
    }

    void simulateSubstep()
    {
        using FEMType = FiniteElementMesh<T, d>;        

        const int nParticles = m_particleX.size();

        // Construct initial guess for next-timestep
        //   Positions -> Same as last timestep
        
        // Overwrite boundary conditions with desired values

        setBoundaryConditions();
        
        // Solve for everything else using Conjugate Gradients

        std::vector<VectorType> dx(nParticles, VectorType::Zero());
        std::vector<VectorType> rhs(nParticles, VectorType::Zero());
        std::vector<VectorType> q(nParticles, VectorType::Zero());
        std::vector<VectorType> s(nParticles, VectorType::Zero());
        std::vector<VectorType> r(nParticles, VectorType::Zero());
        CGVectorWrapper<VectorType> dxWrapper(dx);
        CGVectorWrapper<VectorType> rhsWrapper(rhs);
        CGVectorWrapper<VectorType> qWrapper(q);
        CGVectorWrapper<VectorType> sWrapper(s);
        CGVectorWrapper<VectorType> rWrapper(r);
        CGSystemWrapper<VectorType, FEMType> systemWrapper(*this);
        
        addElasticForce(rhs);
        clearConstrainedParticles(rhs);

        ConjugateGradient<T>::Solve(systemWrapper,
            dxWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper,
            1e-4, 50);

        // Apply corrections to positions and velocities

        for(int p = 0; p < nParticles; p++)
            m_particleX[p] += dx[p];
    }

    void simulateFrame(const int frame)
    {
        m_stepDt = m_frameDt / (T) m_subSteps;

        for(int step = 1; step <= m_subSteps; step++){
            m_stepEndTime = m_frameDt * (T) (frame-1) + m_stepDt * (T) step;
            simulateSubstep();
        }
    }

    virtual void clearConstrainedParticles(std::vector<VectorType>& x) {}
    virtual void setBoundaryConditions() {}
};