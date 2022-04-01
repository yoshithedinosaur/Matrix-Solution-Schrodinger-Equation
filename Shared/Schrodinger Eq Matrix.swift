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
// where H is a square hermitian operator, psi is a column vector, and E is a scalar.
// The solutions to this system can be found by using linear perturbation theory.
// The wave functions are the sums of the coefficients
// The energy solution are the diagonals of the Hamiltonian matrix diagonalized

class SchrodingerEqMat: ObservableObject {
    
    let integrator = Integrator()
    let squareWell = OneDSchrodinger()
    //var ijSquareWellWaveFunctions: [[Double]] = []   // Stores an array of arrays square well solutions (analitically know so we can hard code it in)
    //var nSquareWellEnergy: [Double] = []    // energies of the square well
    //var iSquareWellPotential: [Double] = []
    //var hamiltonianMatrix: [[Double]] = []  // stores a matrix
    //var perturbationTerm: [Double] = []
    
    /// constructHamiltonian:
    /// -constructs the hermitian operator used for finding eigen vectors and eigen values
    /// elements of H: <psi_i| H_ij |psi_j>
    /// the diagonals are the energies of the square well + integrated psi*V*psi
    /// off diagonals just psi*V*psi
    /// parameters:
    /// -potentialType: the potential well type to call
    /// -boxLength: size of potential well
    /// -xStep: step size of potentials
    /// -matrixSize: the size of the hamiltonian, also the highest energy level to calculate En
    func constructHamiltonian(potentialType: String, boxLength: Double, xStep: Double, matrixSize: Int) -> [[Double]] {
        var ijKronecker = 0.0
        var hamiltonianMatrix = [[Double]](repeating: [Double](repeating: 0.0, count: matrixSize), count: matrixSize)
        
        /*
        for i in 0..<matrixSize {
            for j in 0..<matrixSize {
                perturbationTerm[i] = ijSquareWellWaveFunctions[i][j] * iSquareWellPotential[i] * ijSquareWellWaveFunctions[i][j]
            }
        }
        */
        
        //Constructs the Hamiltonian with indices H_ij
        //diagonal elements will have the square well (unperturbed) component + the perturbation component
        for i in 1..<matrixSize+1 {
            for j in 1..<matrixSize+1 {
                if i == j {
                    ijKronecker = 1.0
                } else {
                    ijKronecker = 0.0
                }
                hamiltonianMatrix[i-1][j-1] = ijKronecker * squareWell.squareWellEn(quantNumb: i, boxLength: boxLength) + integrator.perturbationIntegration(potentialType: potentialType, quantumNumbi: i, quantumNumbj: j, boxLength: boxLength, xStep: xStep)
            }
        }
        
        print("\(hamiltonianMatrix)")
        return hamiltonianMatrix
    }
    
    func findWaveFunction(coeffArray: [Double], boxLength: Double, xStep: Double, matrixSize: Int) -> [plotDataType]{
        //let en = squareWell.squareWellEn(quantNumb: quantumNumb, boxLength: boxLength)
        var contentArray = [plotDataType]()
        //var fullPsi: [Double] = []
        var psi: [Double] = []
        //var ek = 0.0
        
        var count: Int = 0
        
        for xVal in stride(from: 0.0, to: boxLength, by: xStep) {
            psi.append(0.0)
            
            for k in 1..<matrixSize+1 {
                    //ek = squareWell.squareWellEn(quantNumb: k, boxLength: boxLength)
                    
                    psi[count] += coeffArray[k-1] * squareWell.squareWellPsi(boxLength: boxLength, quantumNumb: k, xPosition: xVal)
            }
            /*
            psi[count] += coeffArray[0] * squareWell.squareWellPsi(boxLength: boxLength, quantumNumb: 1, xPosition: xVal)
            psi[count] += coeffArray[1] * squareWell.squareWellPsi(boxLength: boxLength, quantumNumb: 2, xPosition: xVal)
            */
            //fullPsi.append(squareWell.squareWellPsi(boxLength: boxLength, quantumNumb: quantumNumb, xPosition: xVal) + correctionPsi[Int(xVal/xStep)])
            
            contentArray.append([.X: xVal, .Y: psi[count]])
            count += 1
        }
        
        
        return contentArray
    }
    
    func findEigenEnergies(realStartingArray: [[Double]]) -> [(eigenEnergy: Double, coeffArray: [Double])] {
        let N = Int32(realStartingArray.count)
        
        let flatArray :[Double] = pack2dArray(arr: realStartingArray, rows: Int(N), cols: Int(N))
        
        let eigenArray = calculateEigenvalues(arrayForDiagonalization: flatArray)
        let sortedArray = eigenArray.sorted(by: {$0.eigenEnergy < $1.eigenEnergy})
        
        //print("\(sortedArray)")
        
        return sortedArray
    }
    
    
    
    /// findEigens
    /// -finds diagonalizes the hamiltonian matrix using dgeev_
    /// calculateEigenvalues
    ///
    /// - Parameter arrayForDiagonalization: linear Column Major FORTRAN Array for Diagonalization
    /// - Returns: String consisting of the Eigenvalues and Eigenvectors
    func calculateEigenvalues(arrayForDiagonalization: [Double]) -> [(eigenEnergy: Double, coeffArray: [Double])] {
        /* Integers sent to the FORTRAN routines must be type Int32 instead of Int */
        //var N = Int32(sqrt(Double(startingArray.count)))
        
        var returnEnergyArray: [Double] = []
        
        
        var N = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N2 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N3 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N4 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        
        var flatArray = arrayForDiagonalization
        
        var error : Int32 = 0
        var lwork = Int32(-1)
        // Real parts of eigenvalues
        var wr = [Double](repeating: 0.0, count: Int(N))
        // Imaginary parts of eigenvalues
        var wi = [Double](repeating: 0.0, count: Int(N))
        // Left eigenvectors
        var vl = [Double](repeating: 0.0, count: Int(N*N))
        // Right eigenvectors
        var vr = [Double](repeating: 0.0, count: Int(N*N))
        
        
        /* Eigenvalue Calculation Uses dgeev */
        /*   int dgeev_(char *jobvl, char *jobvr, Int32 *n, Double * a, Int32 *lda, Double *wr, Double *wi, Double *vl,
         Int32 *ldvl, Double *vr, Int32 *ldvr, Double *work, Int32 *lwork, Int32 *info);*/
        
        /* dgeev_(&calculateLeftEigenvectors, &calculateRightEigenvectors, &c1, AT, &c1, WR, WI, VL, &dummySize, VR, &c2, LWork, &lworkSize, &ok)    */
        /* parameters in the order as they appear in the function call: */
        /* order of matrix A, number of right hand sides (b), matrix A, */
        /* leading dimension of A, array records pivoting, */
        /* result vector b on entry, x on exit, leading dimension of b */
        /* return value =0 for success*/
        
        
        
        /* Calculate size of workspace needed for the calculation */
        
        var workspaceQuery: Double = 0.0
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspaceQuery, &lwork, &error)
        
        print("Workspace Query \(workspaceQuery)")
        
        /* size workspace per the results of the query */
        
        var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
        lwork = Int32(workspaceQuery)
        
        /* Calculate the size of the workspace */
        
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspace, &lwork, &error)
        
        var eigenSystem: (energy: Double, coeff: [Double]) = (energy: 0.0, coeff: [0.0])
        var returningEigenSystems: [(energy: Double, coeff: [Double])] = []
        
        if (error == 0)
        {
            
            
            
            for index in 0..<wi.count      /* transform the returned matrices to eigenvalues and eigenvectors */
            
            {
                if (wi[index]>=0.0)
                {
                    returnEnergyArray.append(wr[index])
                }
                else
                {
                    returnEnergyArray.append(wr[index])
                }
                
                /* To Save Memory dgeev returns a packed array if complex */
                /* Must Unpack Properly to Get Correct Result
                 
                 VR is DOUBLE PRECISION array, dimension (LDVR,N)
                 If JOBVR = 'V', the right eigenvectors v(j) are stored one
                 after another in the columns of VR, in the same order
                 as their eigenvalues.
                 If JOBVR = 'N', VR is not referenced.
                 If the j-th eigenvalue is real, then v(j) = VR(:,j),
                 the j-th column of VR.
                 If the j-th and (j+1)-st eigenvalues form a complex
                 conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
                 v(j+1) = VR(:,j) - i*VR(:,j+1). */
                
                //print("\(vr)")
                
                var returnCoeffArray: [Double] = []
                
                for j in 0..<N
                {
                    //if(wi[index]==0)
                    //{
                        
                    let bob = Int(index)*(Int(N))+Int(j)
                        returnCoeffArray.append(vr[bob])
                    if bob == 10 {
                        print("\(vr[bob])")
                        print("hello")
                    }
                    
                }
                
                eigenSystem = (energy: returnEnergyArray[index], coeff: returnCoeffArray)
                
                returningEigenSystems.append(eigenSystem)
                
                print("\(returningEigenSystems)")
            }
        }
        
    //    let unpackedReturnCoeffArray = unpack2dArray(arr: returnCoeffArray, rows: Int(N), cols: Int(N))
        var returnArray: [(eigenEnergy: Double, coeffArray: [Double])] = []
        for i in 0..<Int(N) {
            returnArray.append((eigenEnergy: returnEnergyArray[i], coeffArray: returningEigenSystems[i].coeff))
        }
        
        return returnArray
    }
    
    /// pack2DArray
    /// Converts a 2D array into a linear array in FORTRAN Column Major Format
    ///
    /// - Parameters:
    ///   - arr: 2D array
    ///   - rows: Number of Rows
    ///   - cols: Number of Columns
    /// - Returns: Column Major Linear Array
    func pack2dArray(arr: [[Double]], rows: Int, cols: Int) -> [Double] {
        var resultArray = Array(repeating: 0.0, count: rows*cols)
        for Iy in 0...cols-1 {
            for Ix in 0...rows-1 {
                let index = Iy * rows + Ix
                resultArray[index] = arr[Ix][Iy]
            }
        }
        return resultArray
    }
    
    /// unpack2DArray
    /// Converts a linear array in FORTRAN Column Major Format to a 2D array in Row Major Format
    ///
    /// - Parameters:
    ///   - arr: Column Major Linear Array
    ///   - rows: Number of Rows
    ///   - cols: Number of Columns
    /// - Returns: 2D array
    func unpack2dArray(arr: [Double], rows: Int, cols: Int) -> [[Double]] {
        var resultArray = [[Double]](repeating:[Double](repeating:0.0 ,count:rows), count:cols)
        for Iy in 0...cols-1 {
            for Ix in 0...rows-1 {
                let index = Iy * rows + Ix
                resultArray[Ix][Iy] = arr[index]
            }
        }
        return resultArray
    }
    
    
    
}

