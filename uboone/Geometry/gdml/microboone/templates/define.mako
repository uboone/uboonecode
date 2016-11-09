<define>
   <rotation name="rPlus30AboutX"  unit="deg" x="30"  y="0"   z="0"/>
   <rotation name="rPlus60AboutX"  unit="deg" x="60"  y="0"   z="0"/>
   <rotation name="rPlus90AboutX"  unit="deg" x="90"  y="0"   z="0"/>
   <rotation name="rMinus90AboutX"  unit="deg" x="-90"  y="0"   z="0"/>
   <rotation name="rPlusUVAngleAboutX"  unit="deg" x="${90.+ attributes['UVAngle']}" y="0"   z="0"/>
   <rotation name="rPlus150AboutX"      unit="deg" x="150" y="0"   z="0"/>
   <rotation name="rPlus180AboutX"      unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rMinusUVAngleAboutX" unit="deg" x="-30" y="0"   z="0"/>
   <rotation name="rPlus30AboutY"  unit="deg" x="0"   y="30"  z="0"/>
   <rotation name="rPlus60AboutY"  unit="deg" x="0"   y="60"  z="0"/>
   <rotation name="rPlus90AboutY"  unit="deg" x="0"   y="90"  z="0"/>
   <rotation name="rPlus180AboutY" unit="deg" x="0"   y="180" z="0"/>
   <rotation name="rMinus90AboutY" unit="deg" x="0"   y="-90" z="0"/>
   <rotation name="rPlus90AboutZ"  unit="deg" x="0"   y="0"   z="90"/>
   <rotation name="rMinus90AboutZ"  unit="deg" x="0"   y="0"   z="-90"/>
   <rotation name="rPlus180AboutZ"      unit="deg" x="0"   y="0"   z="180"/>
   <rotation name="rMinus180AboutZ"     unit="deg" x="0"   y="0"   z="-180"/>
   <rotation name="rMinus90AboutYPlus180AboutZ" unit="deg" x="0" y="-90" z="180"/>
   <rotation name="rMinus90AboutYMinus90AboutZ" unit="deg" x="0" y="-90" z="-90"/>
   <rotation name="rPlus90AboutYPlus180AboutZ" unit="deg" x="0" y="90" z="180"/>
   <rotation name="rMinus90AboutYPlus90AboutZ" unit="deg" x="0" y="-90" z="90"/>
   <rotation name="rPlus90AboutYMinus90AboutZ" unit="deg" x="0" y="90" z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutZ"  unit="deg" x="90" y="0"   z="90"/>
   <rotation name="rPlus90AboutXPlus180AboutZ" unit="deg" x="90" y="0"   z="180"/>
   <rotation name="rPlus90AboutXMinus90AboutY" unit="deg" x="90" y="-90" z="0"/>
   <rotation name="rPlus90AboutXMinus90AboutZ" unit="deg" x="90" y="0"   z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutY"  unit="deg"  x="90" y="90" z="0"/>
   <rotation name="rPMTRotation1"  unit="deg" x="90"  y="270"   z="0"/>
   <rotation name="r39degFix" unit="deg" x="0" y="39" z="0"/>
   <position name="posCenter" unit="mm" x="0" y="0" z="0"/>

   <position name="posEndCap1" unit="cm" x="0" y="0" z="${attributes['EndcapZcenter']}"/>
   <position name="posEndCap2" unit="cm" x="0" y="0" z="-${attributes['EndcapZcenter']}"/>
</define>
