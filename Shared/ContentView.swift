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
    var body: some View {
        
        HStack {
            
            VStack {
                Text("Potential Selector")
                
                Text("Eigen energies")
            }
            
            VStack {
                
            }
            
        }
        
    }
    
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
