function testMutantChristoffel()
c = makeTetragonalMaterialTensor(5.6e10, 5.145e10, 2.2e10, 10.6e10, 6.6e10, 2.65e10);
[mutantPoly, mutantMatrix] = constructEquation(c, 5.96e3);
mutantChristoffel(mutantPoly, mutantMatrix, 0.00001, 0.00001)