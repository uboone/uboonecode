<box name="DetEnclosureOLDbig" lunit="cm" x="${attributes['DetEnclosureWidth']+1}" y="${attributes['DetEnclosureHeight']+1}" z="${attributes['DetEnclosureLength']+1}" />
<box name="DetEnclosureOLDsmall" lunit="cm" x="${attributes['DetEnclosureWidth']}" y="${attributes['DetEnclosureHeight']}" z="${attributes['DetEnclosureLength']}" />
<tube name="DetEnclosureNEW" lunit="cm" rmax="${292*2.54-0.1}" z="${attributes['DetEnclosureHeight']+10}" aunit="deg" deltaphi="360" />
<tube name="DetEnclosureAboveGrade" rmax="${292*2.54-.5}" z="${(29*12+5)*2.54-83*2}" deltaphi="360" lunit="cm" aunit="deg"/>

 <subtraction name="DetSub0">
 <first ref="DetEnclosureOLDbig"/> <second ref="DetEnclosureNEW"/>
 <position name="posDetUnion0" unit="cm" x="0" y="0" z="0"/>
 <rotation name="rotDetUnion0" unit="deg" x="90" y="0" z="0"/>
 </subtraction>

 <subtraction name="DetEnclosure0">
 <first ref="DetEnclosureOLDsmall"/> <second ref="DetSub0"/>
 <position name="posDetEnclosure2" unit="cm" x="0" y="0" z="0"/>
 </subtraction>

 <union name="DetEnclosure">
 <first ref="DetEnclosure0"/> <second ref="DetEnclosureAboveGrade"/>
 <position name="posDetEnclosure3" unit="cm" x="0" y="${attributes['DetEnclosureHeight']/2+(29*12+5)*2.54/2-83}" z="0"/>
 <rotation name="rotDetUnion1" unit="deg" x="90" y="0" z="0"/>
 </union>
