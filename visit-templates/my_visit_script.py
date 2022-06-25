DefineScalarExpression("FeCoNiX","fe56+co56+ni56+x56")
DefineScalarExpression("vr_km","vex/1.e8")

sph2car = TransformAttributes()
sph2car.transformType = sph2car.Coordinate
sph2car.inputCoordSys = sph2car.Spherical
sph2car.outputCoordSys = sph2car.Cartesian
sph2car.vectorTransformMethod = sph2car.AsDirection
sph2car.transformVectors = 1

iso = IsosurfaceAttributes()
iso.contourValue = (0.03)
iso.contourMethod = iso.Value
iso.variable = "FeCoNiX"

AddPlot("Pseudocolor","vr_km")
vr = PseudocolorAttributes()
vr.scaling = vr.Linear
vr.limitsMode = vr.CurrentPlot
vr.centering = vr.Nodal
vr.colorTableName = "RdYlBu"
vr.invertColorTable = 1
SetPlotOptions(vr)
AddOperator("Transform")
SetOperatorOptions(sph2car)
AddOperator("Isosurface")
SetOperatorOptions(iso)

################################################################

DefineScalarExpression("FeCoNiX","fe56+co56+ni56+x56")
DefineScalarExpression("den FeCoNiX","den*(fe56+co56+ni56+x56)")
DefineScalarExpression("den Ti","den*ti44")
DefineScalarExpression("den H","den*p")
DefineScalarExpression("den He","den*he4")
DefineScalarExpression("vr_km","vex/1.e8")
sph2car = TransformAttributes()
sph2car.transformType = sph2car.Coordinate
sph2car.inputCoordSys = sph2car.Spherical
sph2car.outputCoordSys = sph2car.Cartesian
sph2car.vectorTransformMethod = sph2car.AsDirection
sph2car.transformVectors = 1

iso = IsosurfaceAttributes()
iso.contourValue = (2e-16)
iso.contourMethod = iso.Value
iso.variable = "default"
AddPlot("Pseudocolor","den Ti")
vr = PseudocolorAttributes()
vr.scaling = vr.Linear
vr.limitsMode = vr.CurrentPlot
vr.centering = vr.Nodal
vr.colorTableName = "Blues"
vr.invertColorTable = 1
SetPlotOptions(vr)
AddOperator("Transform")
SetOperatorOptions(sph2car)
AddOperator("Isosurface")
SetOperatorOptions(iso)

iso = IsosurfaceAttributes()
iso.contourValue = (0.03)
iso.contourMethod = iso.Value
iso.variable = "FeCoNiX"
AddPlot("Pseudocolor","vr_km")
vr = PseudocolorAttributes()
vr.scaling = vr.Linear
vr.limitsMode = vr.CurrentPlot
vr.centering = vr.Nodal
vr.colorTableName = "RdYlBu"
vr.invertColorTable = 1
SetPlotOptions(vr)
AddOperator("Transform")
SetOperatorOptions(sph2car)
AddOperator("Isosurface")
SetOperatorOptions(iso)

iso = IsosurfaceAttributes()
iso.contourValue = (3e-14)
iso.contourMethod = iso.Value
iso.variable = "default"
AddPlot("Pseudocolor","den He")
vr = PseudocolorAttributes()
vr.scaling = vr.Linear
vr.limitsMode = vr.CurrentPlot
vr.centering = vr.Nodal
vr.colorTableName = "Greens"
vr.invertColorTable = 1
SetPlotOptions(vr)
AddOperator("Transform")
SetOperatorOptions(sph2car)
AddOperator("Isosurface")
SetOperatorOptions(iso)

iso = IsosurfaceAttributes()
iso.contourValue = (1e-15)
iso.contourMethod = iso.Value
iso.variable = "default"
AddPlot("Pseudocolor","den He")
vr = PseudocolorAttributes()
vr.scaling = vr.Linear
vr.limitsMode = vr.CurrentPlot
vr.centering = vr.Nodal
vr.colorTableName = "Reds"
vr.invertColorTable = 1
SetPlotOptions(vr)
AddOperator("Transform")
SetOperatorOptions(sph2car)
AddOperator("Isosurface")
SetOperatorOptions(iso)

iso = IsosurfaceAttributes()
iso.contourValue = (1e-18)
iso.contourMethod = iso.Value
iso.variable = "default"
AddPlot("Pseudocolor","den H")
vr = PseudocolorAttributes()
vr.scaling = vr.Linear
vr.limitsMode = vr.CurrentPlot
vr.centering = vr.Nodal
vr.colorTableName = "Oranges"
vr.invertColorTable = 1
SetPlotOptions(vr)
AddOperator("Transform")
SetOperatorOptions(sph2car)
AddOperator("Isosurface")
SetOperatorOptions(iso)

ResetView()
DrawPlots()
