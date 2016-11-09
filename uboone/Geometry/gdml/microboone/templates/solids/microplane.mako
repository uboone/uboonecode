% for i in range(attributes['NumberWiresPerEdge']):
<tube name="TPCWire${i}"
  rmax="${0.5*attributes['TPCWireThickness']}"
  z="${attributes['TPCYWirePitch']*(i+1.) * 2. -0.1}"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
% endfor
<tube name="TPCWireCommon"
  rmax="${0.5*attributes['TPCWireThickness']}"
  z="${attributes['TPCWirePlaneWidth']/attributes['CosUVAngle'] - 0.1}"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TPCPlane"
  x="${attributes['TPCWirePlaneThickness']}"
  y="${attributes['TPCWirePlaneWidth']}"
  z="${attributes['TPCWirePlaneLength']+0.35}"
  lunit="cm"/>
