//
//  Wave Functions.swift
//  Matrix Solution Schrodinger Equation
//
//  Created by Yoshinobu Fujikake on 3/11/22.
//

import Foundation

// Last assignment use to get square well wave functions psi_i
class OneDSchrodinger: ObservableObject {
    
    @Published var oneDPotentialArray: [(xCoord: Double, Potential: Double)] = []
    @Published var oneDPotentialXArray: [Double] = []
    @Published var oneDPotentialYArray: [Double] = []
    @Published var energyVal = 0.0
    
    let hbar2overm = 7.63 //In units of eV * Å^2
    
    ///squareWellWaveFunction: known solution to the infinite square well potential
    //                  _
    //                 /2     /n pi x\
    //  psi (x)  =  | / - sin |------|
    //     n        |/  L     \   L  /
    /// parameters:
    /// -boxLength: the size of the square well
    /// -quantumNumb: quantum number of the wave function > 0
    /// -xPosition: the position in the box
    func squareWellPsi(boxLength: Double, quantumNumb: Int, xPosition: Double) -> Double {
        var psiX: Double
        
        psiX = sqrt(2.0/boxLength) * sin(Double(quantumNumb) * Double.pi * xPosition / boxLength)

        
        return psiX
    }
    
    
    ///squareWellEn: gives the En energy of the square well for quantum number n
    //           2  2   2
    //      hbar   n  pi
    // E  = -------------
    //  n           2
    //         2 m L
    ///parameters:
    ///-quantNumb: quantum number n
    ///-boxLength: size of square well box
    func squareWellEn(quantNumb: Int, boxLength: Double) -> Double {
        var energyN: Double
        
        energyN = hbar2overm/2.0 * pow(Double(quantNumb), 2.0) * pow(Double.pi, 2.0) / pow(boxLength, 2.0)
        
        return energyN
    }
    
    
}
