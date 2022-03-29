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
    
    func perturbationIntegration(potentialType: String, quantumNumbi: Int, quantumNumbj: Int, boxLength: Double, xStep: Double) -> Double {
        var integration = 0.0
        var pointsCount: Int = 0
        potentials.getPotential(potentialType: potentialType, xMin: 0.0, xMax: boxLength, xStep: xStep)
        
        for xVal in stride(from: 0.0, through: boxLength, by: xStep) {
            integration += oneDSchrodinger.squareWellPsi(boxLength: boxLength, quantumNumb: quantumNumbi, xPosition: xVal) * oneDSchrodinger.oneDPotentialYArray[Int(xVal/xStep)] * oneDSchrodinger.squareWellPsi(boxLength: boxLength, quantumNumb: quantumNumbj, xPosition: xVal)
            pointsCount += 1
        }
        
        integration = integration/Double(pointsCount)
        
        return integration
    }
    
}
