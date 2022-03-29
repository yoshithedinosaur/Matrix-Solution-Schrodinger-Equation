//
//  Schrodinger Eq Matrix.swift
//  Matrix Solution Schrodinger Equation
//
//  Created by Yoshinobu Fujikake on 3/11/22.
//

import Foundation
import Accelerate

// ^
// H psi = E psi
// where H is a square hermitian operator, psi is a column vector, and E is a scalar
// H maybe can be found via jacobian since H = hbar2overm/2 partialx^2 + V

// Solution: 0 = det(H - lambda * I)
// Find lambda
// lambda finds psi

// Matrix manipulations we have:
// dgeev_
// -query once to get workspace
// -query again with calculated workspace and return a bunch of values including, real and imaginary eigenvalues, and eigenvectors
// Note: wavefunction is the sum of coefficients * eigenvectors

class SchrodingerEqMat: ObservableObject {
    
    let integrator = Integrator()
    let squareWell = OneDSchrodinger()
    //var ijSquareWellWaveFunctions: [[Double]] = []   // Stores an array of arrays square well solutions (analitically know so we can hard code it in)
    //var nSquareWellEnergy: [Double] = []    // energies of the square well
    //var iSquareWellPotential: [Double] = []
    var hamiltonianMatrix: [[Double]] = []  // stores a matrix
    //var perturbationTerm: [Double] = []
    
    /// constructHamiltonian
    /// -constructs the hermitian operator used for finding eigen vectors and eigen values
    /// elements of H: <psi_i| H_ij |psi_j>
    /// the diagonals are the energies of the square well + integrated psi*V*psi
    /// off diagonals just psi*V*psi
    func constructHamiltonian(potentialType: String, boxLength: Double, xStep: Double, matrixSize: Int) -> [[Double]] {
        var ijKronecker = 0.0
        
        /*
        for i in 0..<matrixSize {
            for j in 0..<matrixSize {
                perturbationTerm[i] = ijSquareWellWaveFunctions[i][j] * iSquareWellPotential[i] * ijSquareWellWaveFunctions[i][j]
            }
        }
        */
        
        //Constructs the Hamiltonian with indices H_ij
        //diagonal elements will have the square well (unperturbed) component + the perturbation component
        for i in 0..<matrixSize {
            for j in 0..<matrixSize {
                if i == j {
                    ijKronecker = 1.0
                } else {
                    ijKronecker = 0.0
                }
                hamiltonianMatrix[i][j] = ijKronecker * squareWell.squareWellEn(quantNumb: i, boxLength: boxLength) + integrator.perturbationIntegration(potentialType: potentialType, quantumNumbi: i, quantumNumbj: j, boxLength: boxLength, xStep: xStep)
            }
        }
        
        return hamiltonianMatrix
    }
    
    
    
    /// findEigens
    /// -finds all eigen values and eigen vectors using dgeev_
    func findEigens() {
        
    }
    
    
    
}

