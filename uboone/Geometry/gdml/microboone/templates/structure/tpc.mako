% if attributes['cryostat_on']:
  <volume name="volTPCActive">
   <materialref ref="LAr"/>
   <solidref ref="TPCActive"/>
  </volume>
  <volume name="volTPC">
   <materialref ref="LAr"/>
   <solidref ref="TPC"/>
   <physvol>
    <volumeref ref="volCathodePlate"/>
    <position name="posCathodePlate" unit="cm"
    x="${0.5*(attributes['TPCDepth']-attributes['CathodeWidth'])}"
     y="0" z="0"/>
   </physvol>
   <physvol>
    <volumeref ref="volTPCPlaneVert"/>
    <position name="posTPCPlaneVert" unit="cm"
    x="${-0.5*attributes['TPCActiveDepth']-2*attributes['TPCWirePlaneSpacing']}"
    y="0" z="0"/>
   </physvol>
   <physvol>
    <volumeref ref="volTPCPlane"/>
    <position name="posTPCPlane" unit="cm"
    x="${-0.5*attributes['TPCActiveDepth']-attributes['TPCWirePlaneSpacing']}"
    y="0" z="0"/>
   </physvol>
   <physvol>
    <volumeref ref="volTPCPlane"/>
    <position name="posTPCPlane2" unit="cm"
    x="${-0.5*attributes['TPCActiveDepth']}"
    y="0" z="0"/>
    <rotationref ref="rPlus180AboutY"/>
   </physvol>
   % for space in range(0,attributes['field_cage_width']/2, 4*attributes['field_cage_loop_interval']):
     <physvol>
      <volumeref ref="volFieldCageTubeTop"/>
      <position name="posFieldCageTubeTopA${loop.index+1}" unit="cm"
        x="${space+2}"
        y="${0.5*(attributes['FieldCageLoopHeight']-2.*attributes['FieldCageTubeRadius'])}"
        z="0"/>
     </physvol>
     <physvol>
      <volumeref ref="volFieldCageTubeBot"/>
      <position name="posFieldCageTubeBotA${loop.index+1}" unit="cm"
        x="${space+2}"
         y="${-0.5*(attributes['FieldCageLoopHeight']-2.*attributes['FieldCageTubeRadius'])}" z="0"/>
     </physvol>
     <physvol>
      <volumeref ref="volFieldCageTubeFront"/>
      <position name="posFieldCageTubeFrontA${loop.index+1}" unit="cm"
        x="${space+2}"
         y="0" z="${0.5*(attributes['FieldCageLoopLength']-2.*attributes['FieldCageTubeRadius'])}"/>
      <rotation name="rFieldCageVertPlusA${loop.index+1}" unit="deg" x="90" y="0" z="0"/>
     </physvol>
     <physvol>
      <volumeref ref="volFieldCageTubeBack"/>
      <position name="posFieldCageTubeBackA${loop.index+1}" unit="cm"
        x="${space+2}"
        y="0" z="${-0.5*(attributes['FieldCageLoopLength']-2.*attributes['FieldCageTubeRadius'])}"/>
      <rotation name="rFieldCageVertMinusA${loop.index+1}" unit="deg" x="-90" y="0" z="0"/>
     </physvol>
     <physvol>
      <volumeref ref="volFieldCageTubeTop"/>
      <position name="posFieldCageTubeTopB${loop.index+1}" unit="cm"
        x="-${space+2}"
        y="${0.5*(attributes['FieldCageLoopHeight']-2.*attributes['FieldCageTubeRadius'])}" z="0"/>
     </physvol>
     <physvol>
      <volumeref ref="volFieldCageTubeBot"/>
      <position name="posFieldCageTubeBotB${loop.index+1}" unit="cm"
        x="-${space+2}"
        y="${-0.5*(attributes['FieldCageLoopHeight']-2.*attributes['FieldCageTubeRadius'])}" z="0"/>
     </physvol>
     <physvol>
      <volumeref ref="volFieldCageTubeFront"/>
      <position name="posFieldCageTubeFrontB${loop.index+1}" unit="cm"
        x="-${space+2}"
        y="0" z="${0.5*(attributes['FieldCageLoopLength']-2.*attributes['FieldCageTubeRadius'])}"/>
      <rotation name="rFieldCageVertFrontB${loop.index+1}" unit="deg" x="90" y="0" z="0"/>
     </physvol>
     <physvol>
      <volumeref ref="volFieldCageTubeBack"/>
      <position name="posFieldCageTubeBackB${loop.index+1}" unit="cm"
        x="-${space+2}"
        y="0" z="${-0.5*(attributes['FieldCageLoopLength']-2.*attributes['FieldCageTubeRadius'])}"/>
      <rotation name="rFieldCageVertBackB${loop.index+1}" unit="deg" x="-90" y="0" z="0"/>
     </physvol>
   %endfor
     <physvol>
      <volumeref ref="volTPCActive"/>
      <position name="posTPCActive" unit="cm" x="-${attributes['aTPC_xoffset']}" y="0" z="0"/>
     </physvol>
  </volume>
%endif
