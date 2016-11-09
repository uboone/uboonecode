%if attributes['cryostat_on']:

<volume name="volUllage">
  <materialref ref="Argon_gas_87K"/>
  <solidref ref="Ullage"/>
</volume>
<volume name="volSteelVessel">
 <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
 <solidref ref="SteelVessel"/>
</volume>
<volume name="volaTopBeam">
   <materialref ref="G10"/>
   <solidref ref="aTopBeam"/>
</volume>
<volume name="volaTopCross">
   <materialref ref="G10"/>
   <solidref ref="aTopCross"/>
</volume>
<volume name="volaTopCrossOuter">
   <materialref ref="G10"/>
   <solidref ref="aTopCrossOuter"/>
</volume>
<volume name="volaSideCross">
   <materialref ref="G10"/>
   <solidref ref="aSideCross"/>
</volume>
<volume name="volFrame">
 <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
 <solidref ref="Frame4"/>
</volume>
<volume name="volCrossBeam">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CrossBeam"/>
</volume>

<volume name="volCryostat">
 <materialref ref="LAr"/>
 <solidref ref="Cryostat"/>
 <physvol>
  <volumeref ref="volUllage"/>
 </physvol>
 <physvol>
  <volumeref ref="volSteelVessel"/>
  <position name="posSteelVessel" unit="cm" x="0" y="0" z="0"/>
 </physvol>
<physvol>
   <volumeref ref="volFrame"/>
   <position name="posFrame2" unit="cm" x="${-256./2.-8}" y="0" z="0"/>
 </physvol>
 <physvol>
     <volumeref ref="volaSideCross"/>
     <position name="posaSideCross1" unit="cm" x="0" y="-45.63" z="${1050.19/2.+5}"/>
     <rotationref ref="rPlus90AboutX"/>
 </physvol>
 <physvol>
     <volumeref ref="volaSideCross"/>
     <position name="posaSideCross2" unit="cm" x="0" y="-45.63" z="${-1050.19/2.-5}"/>
     <rotationref ref="rPlus90AboutX"/>
 </physvol>
 %for i in range(3):
   <physvol>
       <volumeref ref="volaTopCross"/>
       <position name="posaTopCross${i}" unit="cm" x="0" y="${256./2.+2.5*2.54}" z="${-5./2.*111.64 + 2.*111.64*i}"/>
   </physvol>
   <physvol>
       <volumeref ref="volaTopCross"/>
       <position name="posaBottomCross${i}" unit="cm" x="0" y="${-256./2.-2.5*2.54}" z="${-5./2.*111.64 +2*111.64*i}"/>
   </physvol>
 %endfor
 %for i in range(2):
   <physvol>
       <volumeref ref="volaTopCrossOuter"/>
       <position name="posaTopCrossOuter${i}" unit="cm" x="0" y="${256./2.+2.5*2.54}" z="${-4*111.64 + 15./2.*111.64*i}"/>
   </physvol>
   <physvol>
       <volumeref ref="volaTopCrossOuter"/>
       <position name="posaBottomCrossOuter${i}" unit="cm" x="0" y="${-256./2.-2.5*2.54}" z="${-4.*111.64 + 15./2.*111.64*i}"/>
   </physvol>
 %endfor
 <physvol>
  <volumeref ref="volTPC"/>
  <position name="posTPC" unit="cm" x="0.0" y="0.97" z="0"/>
 </physvol>
 %if attributes['pmt_on']:
   %for pmt_pos in attributes['pmt_pos']:
   <physvol>
     <volumeref ref="volPMT"/>
     <position name="posPMT${loop.index}" unit="cm" x="${pmt_pos['x']}"
     y="${pmt_pos['y']}" z="${pmt_pos['z']}"/>
     <rotationref ref="rPMTRotation1"/>
   </physvol>
   %endfor
   %for pmt_pos in attributes['pmt_paddle_pos']:
   <physvol>
    <volumeref ref="volPaddle_PMT"/>
    <position name="posPMT${loop.index+32}" unit="cm"
    x="${pmt_pos['x']}"
    y="${pmt_pos['y']}"
    z="${pmt_pos['z']}"/>
   </physvol>
   %endfor
 %endif
</volume>

%endif
