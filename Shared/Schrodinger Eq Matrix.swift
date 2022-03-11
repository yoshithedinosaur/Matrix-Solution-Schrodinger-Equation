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
    var iSquareWellWaveFunctions: [[Double]] = []   // Stores an array of arrays square well solutions (analitically know so we can hard code it in)
    var nSquareWellEnergy: [Double] = []    // energies of the square well
    var iSquareWellPotential: [Double] = []
    var ijHamiltonian: [[Double]] = []  // stores a matrix
    var iPsiVPsi: [Double] = []
    var valuesToIntegrate: [Double] = []
    
    /// constructHamiltonian
    /// -constructs the hermitian operator used for finding eigen vectors and eigen values
    /// elements of H: <psi_i| H_ij |psi_j>
    /// the diagonals are the energies of the square well + integrated psi*V*psi
    /// off diagonals just psi*V*psi
    func constructHamiltonian () {
        var ijKronecker = 0.0
        var psiVPsi = 0.0
        
        for i in 0..<nSquareWellEnergy.count {
            valuesToIntegrate[i] = iPsiVPsi[i] * iSquareWellPotential[i] * iPsiVPsi[i]
        }
        
        for i in 0..<nSquareWellEnergy.count {
            for j in 0..<nSquareWellEnergy.count {
                if i == j {
                    ijKronecker = 1.0
                } else {
                    ijKronecker = 0.0
                }
                psiVPsi = integrator.avgValIntegration(valuesToIntegrate: valuesToIntegrate)
                ijHamiltonian[i][j] = ijKronecker * nSquareWellEnergy[i] + psiVPsi
            }
        }
    }
    
    /// findEigens
    /// -finds all eigen values and eigen vectors using dgeev_
    func findEigens() {
        
    }
    
    
    
}

