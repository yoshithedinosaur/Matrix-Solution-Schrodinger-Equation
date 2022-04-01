//
//  ContentView.swift
//  Shared
//
//  Created by Yoshinobu Fujikake on 3/11/22.
//

import SwiftUI
import CorePlot

typealias plotDataType = [CPTScatterPlotField : Double]

struct ContentView: View {
    @ObservedObject var plotDataModel = PlotDataClass(fromLine: true)
    @ObservedObject private var potentialPlotter = PlotPotentials()
    @ObservedObject private var waveFunctionPlotter = PlotWaveFunctions()
    @ObservedObject private var potentialWells = PotentialWells()
    @ObservedObject private var oneDSchrodinger = OneDSchrodinger()
    @ObservedObject private var schrodingerEqMat = SchrodingerEqMat()
    //@ObservedObject private var energies = EigenEnergyFinder()
    
    @State private var potentialSelect = "Square Well"
    @State private var viewSelect = "Potential"
    @State var energyLevelView = ""
    
    @State var matrixSizeString = "4"
    
    //Choices set up
    let potentialsChoices = ["Square Well", "Linear Well", "Parabolic Well", "Square + Linear Well", "Square Barrier", "Triangle Barrier", "Coupled Parabolic Well", "Coupled Square Well + Field"/*, "Harmonic Oscillator"*/, "Kronig - Penney"/*, "Variable Kronig - Penney", "KP2-a"*/]
    let viewChoices = ["Potential", "Wave Function"]
    @State var eigenEnergyList: [String] = []
    
    
    var body: some View {
        
        HStack {
            
            VStack {
                
                Text("Matrix Size")
                TextField("", text: $matrixSizeString)

                Text("Choose Potential Well")
                Picker("", selection: $potentialSelect) {
                    ForEach(potentialsChoices, id: \.self) {
                        Text($0)
                    }
                }
                
                Text("Choose View")
                Picker("", selection: $viewSelect) {
                    ForEach(viewChoices, id: \.self) {
                        Text($0)
                    }
                }
                /*.onReceive([self.viewSelect].publisher.first()) { potentialSelect in
                    plotPotentialWells(potentialType: potentialSelect)
                }*/ //Figure out later
                
                
                Text("Energy Levels")
                Picker("", selection: $energyLevelView) {
                    ForEach(eigenEnergyList, id: \.self) {
                        Text($0)
                    }
                }
                
                Text("Select potential first, then click 'Print Energies' to get the energy levels to plot the wave functions of. (the wave functions are not correct btw)")
                    .frame(width: 400, height: 50)
            }
            .frame(width: 400, height: 40)
            .padding()
            .frame(width: 400, height: 40)
            .padding()
            
            VStack{
                CorePlot(dataForPlot: $plotDataModel.plotData, changingPlotParameters: $plotDataModel.changingPlotParameters)
                    .setPlotPadding(left: 10)
                    .setPlotPadding(right: 10)
                    .setPlotPadding(top: 10)
                    .setPlotPadding(bottom: 10)
                    .padding()
                            
                Divider()
                            
                HStack{
                    Button("Plot Selected View", action: {self.plotButton(viewChoice: viewSelect, potentialType: potentialSelect)} )
                        .padding()
                
                    Button("Print Energies", action: {self.printButton(potentialType: potentialSelect, matrixSize: Int(matrixSizeString)!)} )
                        .padding()
                    
                }
                
            }
            
        }
        
    }

    func printButton(potentialType: String, matrixSize: Int) {
        
        eigenEnergyList.removeAll()
        let hamiltonianMatrix = schrodingerEqMat.constructHamiltonian(potentialType: potentialType, boxLength: 10.0, xStep: 0.01, matrixSize: matrixSize)
        
        for i in 0..<matrixSize {
            eigenEnergyList.append("n = " + "\(i+1)" + " (E = \(schrodingerEqMat.findEigenEnergies(realStartingArray: hamiltonianMatrix)[i].eigenEnergy))")
        }
        
    }
    
    func plotButton (viewChoice: String, potentialType: String) {
        switch viewChoice {
        case "Potential":
            plotPotentialWells(potentialType: potentialType)
            
        case "Wave Function":
            if let quantumNumb = eigenEnergyList.firstIndex(of: energyLevelView) {
            plotWaveFunction(potentialType: potentialType, quantumNumb: Int(quantumNumb), matrixSize: Int(matrixSizeString)!)
            } else {
                print("ERROR")
            }
            
        default:
            print("Don't know how you got here but ok.")
        }
    }
    
    func plotPotentialWells(potentialType: String){
                
        //pass the plotDataModel to the potentialPlotter
        potentialPlotter.plotDataModel = self.plotDataModel
        //Calculate the new plotting data and place in the plotDataModel
        potentialPlotter.plotWells(potentialType: potentialType)
        
    }
    
    func plotWaveFunction(potentialType: String, quantumNumb: Int, matrixSize: Int) {
        
        //pass the plotDataModel to the waveFunctionPlotter
        waveFunctionPlotter.plotDataModel = self.plotDataModel
        //Calculate the new plotting data and place in the plotDataModel
        waveFunctionPlotter.plotWaveFunction(potentialType: potentialType, quantumNumb: quantumNumb, matrixSize: matrixSize)
        
    }
    
}


struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
