//
//  Average Value Integration.swift
//  Matrix Solution Schrodinger Equation
//
//  Created by Yoshinobu Fujikake on 3/11/22.
//

import Foundation

class Integrator: ObservableObject {
    let oneDSchrodinger = OneDSchrodinger()
    let potentials = PotentialWells()
    
    ///perturbationIntegration: integrates the perturbation term in the Hamiltonian
    //    _
    //   /
    //   |   psi (x) V(x) psi (x) dx
    //  _/      i            j
    ///parameters:
    ///-potentialType: the type of potential describes the perturbation to the square well potential
    ///-quantumNumbi and quantumNumbj: the quantum numbers of the two wave functions to integrate over
    ///-boxLength: size of potential well
    ///-xStep: step size through the potential well
    func perturbationIntegration(potentialType: String, quantumNumbi: Int, quantumNumbj: Int, boxLength: Double, xStep: Double) -> Double {
        var integration = 0.0
        var pointsCount: Int = 0
        let potentialsArray = potentials.getPotential(potentialType: potentialType, xMin: 0.0, xMax: boxLength, xStep: xStep)
        
        for xVal in stride(from: 0.0, through: boxLength, by: xStep) {
            integration += oneDSchrodinger.squareWellPsi(boxLength: boxLength, quantumNumb: quantumNumbi, xPosition: xVal) * potentialsArray[Int(xVal/xStep)][.Y]! * oneDSchrodinger.squareWellPsi(boxLength: boxLength, quantumNumb: quantumNumbj, xPosition: xVal)
            pointsCount += 1
        }
        
        integration = integration/Double(pointsCount)
        
        return integration
    }
    
}
