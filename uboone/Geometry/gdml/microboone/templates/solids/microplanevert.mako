<tube name="TPCWireVert"
  rmax="${0.5*attributes['TPCWireThickness']}"
  z="${attributes['TPCWirePlaneWidth']}"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TPCPlaneVert"
  x="${attributes['TPCWirePlaneThickness']}"
  y="${attributes['TPCWirePlaneWidth']}"
  z="${attributes['TPCWirePlaneLength']}"
  lunit="cm"/>
