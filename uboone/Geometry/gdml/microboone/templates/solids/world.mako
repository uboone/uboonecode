<box name="World"
  lunit="cm"
  x="${attributes['WorldWidth']}"
  y="${attributes['WorldHeight']}"
  z="${attributes['WorldLength']}"/>
<box name="VacuumBox"
lunit="cm"
x="${attributes['WorldWidth']}"
  y="${53000-1800}"
  z="${attributes['WorldLength']}"/>
<tube name="LArTFVacuumSubtract"
  rmax="${310*2.54}"
  z="${(29*12+6)*2.54}"
  deltaphi="360"
  lunit="cm"
  aunit="deg"/>

<tube name="ConcreteWallAboveGrade0"
 rmin="${292*2.54}"
 rmax="${(310)*2.54}"
 z="${(29*12+6)*2.54}"
 deltaphi="360"
 lunit="cm"
 aunit="deg"/>
<tube name="ConcreteWallAboveGradeSegment"
 rmax="${304*2.54}"
 z="${(4*12+7)*2.54}"
 deltaphi="360"
 lunit="cm"
 aunit="deg"/>
<tube name="ConcreteDiscRoof0"
rmax="${304*2.54-31.5}"
z="${13*2.54}"
deltaphi="360"
lunit="cm"
aunit="deg"/>
<box name="ConcreteSteelBeam0"
  lunit="cm"
  x="${48*2.54}"
  y="${2*220*2.54}"
  z="${(5*12+8)*2.54}"/>
<box name="ConcreteSteelBeamSegment1"
  lunit="cm"
  x="${6*2.54+0.1}"
  y="${2*220*2.54 +0.1}"
  z="${34*2.54+0.1}"/>
<box name="ConcreteSteelBeamSegment2"
  lunit="cm"
  x="${12*2.54+0.1}"
  y="${2*230*2.54+0.1}"
  z="${34*2.54+0.1}"/>

<tube name="RoofPlusBars"
rmax="${304*2.54-31}"
z="${68*2.54 + 0.1}"
deltaphi="360"
lunit="cm"
aunit="deg"/>
<box name="SteelGrating"
  lunit="cm"
x="${(23*12*2.54)}"
  y="${.5*2.54}"
z="${(43*12-2)*2.54}"/>
<box name="RemovableRoof"
lunit="cm"
x="${(23*12*2.54+.1)}"
y="${(42*12+10)*2.54}"
z="${(12.5+12.5)*2.54}"/>
<box name="IBeam22_0"
lunit="cm"
x="${4.03*2.54 -0.1}"
y="${(42*12+10)*2.54-0.1}"
z="${12.5*2.54-0.1}"/>
<box name="IBeam22Segment"
lunit="cm"
x="${(4.03-0.26)*2.54/2}"
y="${(43*12+2)*2.54}"
z="${(12.5-2*0.425)*2.54}"/>
<box name="IBeam35_0"
lunit="cm"
x="${6.56*2.54-0.1}"
y="${(23*12)*2.54-0.1}"
z="${12.5*2.54-0.1}"/>
<box name="IBeam35Segment"
lunit="cm"
x="${(6.56-0.3)*2.54/2}"
y="${(23*12)*2.54}"
z="${(12.5-2*0.52)*2.54}"/>

<tube name="Ground"
  rmin="${312*2.54}"
  rmax="${((50*12)+310)*2.54}"
  z="${41*12*2.54}"
  deltaphi="360"
  lunit="cm"
  aunit="deg"/>
<tube name="GroundBottom"
  rmax="${((50*12)+310)*2.54}"
  z="${50*12*2.54}"
  deltaphi="360"
  lunit="cm"
  aunit="deg"/>
<tube name="ConcreteEnclosure"
  rmin="${292*2.54}"
  rmax="${310*2.54}"
  z="${38*12*2.54}"
  deltaphi="360"
  lunit="cm"
  aunit="deg"/>
<tube name="PolystyreneEnclosure"
  rmin="${310*2.54}"
  rmax="${312*2.54}"
  z="${(38*12+36)*2.54}"
  deltaphi="360"
  lunit="cm"
  aunit="deg"/>
<tube name="ConcreteEnclosureBottom"
  rmin="0"
  rmax="${310*2.54}"
  z="${36*2.54}"
  deltaphi="360"
  lunit="cm"
  aunit="deg"/>
<tube name="PolystyreneEnclosureBottom"
 rmax="${292*2.54}"
 z="${2*2.54}"
 deltaphi="360"
 lunit="cm"
 aunit="deg"/>
<tube name="Overburden"
  rmax="${292*2.54}"
  z="${10*12*2.54}"
  deltaphi="360"
  lunit="cm"
  aunit="deg"/>

<subtraction name="VacuumSpace">
 <first ref="VacuumBox"/> <second ref="LArTFVacuumSubtract"/>
 <rotation name="rotLArTFVac" unit="deg" x="90" y="0" z="0"/>
 <position name="posLArTFVac" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="${-(53000-1800)/2+(29*12+6)*2.54/2}" z="${0.5*attributes['TPCWirePlaneLength']}" />
</subtraction>

<subtraction name="ConcreteWallAboveGrade">
  <first ref="ConcreteWallAboveGrade0"/> <second ref="ConcreteWallAboveGradeSegment"/>
<position name="posRemoval0" unit="cm" x="0" y="0" z="${(122+24+3.5)*2.54}"/>
</subtraction>

<subtraction name="ConcreteDiscRoof">
  <first ref="ConcreteDiscRoof0"/> <second ref="RemovableRoof"/>
<position name="posRemoval1" unit="cm" x="0" y="0" z="0"/>
</subtraction>

<subtraction name="ConcreteSteelBeam1">
  <first ref="ConcreteSteelBeam0"/> <second ref="ConcreteSteelBeamSegment1"/>
<position name="posBeamRemoval0" unit="cm" x="${21*2.54}" y="0" z="${17*2.54}"/>
</subtraction>

<subtraction name="ConcreteSteelBeam">
  <first ref="ConcreteSteelBeam1"/> <second ref="ConcreteSteelBeamSegment2"/>
<position name="posBeamRemoval1" unit="cm" x="${-18*2.54}" y="0" z="${17*2.54}"/>
</subtraction>

<subtraction name="IBeam22_1">
  <first ref="IBeam22_0"/> <second ref="IBeam22Segment"/>
<position name="posBeamRemoval2" unit="cm" x="${(4.03 -1.885)*0.5*2.54}" y="0" z="0"/>
</subtraction>

<subtraction name="IBeam22">
  <first ref="IBeam22_1"/> <second ref="IBeam22Segment"/>
<position name="posBeamRemoval3" unit="cm" x="${-(4.03 -1.885)*0.5*2.54}" y="0" z="0"/>
</subtraction>

<subtraction name="IBeam35_1">
  <first ref="IBeam35_0"/> <second ref="IBeam35Segment"/>
<position name="posBeamRemoval4" unit="cm" x="${(6.56 - 3.13)*0.5*2.54}" y="0" z="0"/>
</subtraction>

<subtraction name="IBeam35">
  <first ref="IBeam35_1"/> <second ref="IBeam35Segment"/>
<position name="posBeamRemoval5" unit="cm" x="${-(6.56-3.13)*0.5*2.54}" y="0" z="0"/>
</subtraction>
