%if attributes['cryostat_on']:
  <tube name="CryostatTube"
   rmax="${attributes['CryostatOuterRadius']}"
   z="${attributes['CryostatLength']}"
   deltaphi="360"
   aunit="deg"
   lunit="cm"/>
  <sphere name="CryostatEnd" rmin="0" rmax="${attributes['EndcapRmax']}" deltaphi="360" deltatheta="${attributes['EndcapThetaDeg']}" aunit="deg" lunit="cm"/>
   <union name="CryostatUnion1">
    <first ref="CryostatTube"/>
    <second ref="CryostatEnd"/>
    <positionref ref="posEndCap1"/>
   </union>
   <union name="Cryostat">
    <first ref="CryostatUnion1"/>
    <second ref="CryostatEnd"/>
    <positionref ref="posEndCap2"/>
    <rotationref ref="rPlus180AboutY"/>
   </union>
   <tube name="UllageTube"
    rmax="${attributes['CryostatInnerRadius']}"
    z="${attributes['UllageLength']}"
    deltaphi="360"
    aunit="deg"
    lunit="cm"/>
   <sphere name="UllageEnd" rmin="0" rmax="${attributes['EndcapRmin']}" deltaphi="360" deltatheta="${attributes['UllageThetaDeg']}" aunit="deg" lunit="cm"/>
    <union name="UllageUnion1">
     <first ref="UllageTube"/>
     <second ref="UllageEnd"/>
     <positionref ref="posEndCap1"/>
    </union>
    <union name="UllageUnion2">
     <first ref="UllageUnion1"/>
     <second ref="UllageEnd"/>
     <positionref ref="posEndCap2"/>
     <rotationref ref="rPlus180AboutY"/>
    </union>
    <box name="UllageBox" lunit="cm" x="${2*attributes['CryostatInnerRadius']}" y="${attributes['UllageLevelFromTop']}" z="${attributes['CryostatLength']+attributes['CryostatEndcapLength']}"/>
    <intersection name="Ullage">
     <first ref="UllageUnion2"/>
     <second ref="UllageBox"/>
     <position name="posUllageBox" unit="cm" x="0" y="${attributes['CryostatInnerRadius']-0.5*attributes['UllageLevelFromTop']}" z="0"/>
    </intersection>
    <tube name="SteelTube"
      rmin="${attributes['CryostatInnerRadius']}"
      rmax="${attributes['CryostatOuterRadius']-0.0001}"
      z="${attributes['CryostatLength']}"
      deltaphi="360"
      aunit="deg"
      lunit="cm"/>

    <sphere name="EndCap" rmin="${attributes['EndcapRmin']}" rmax="${attributes['EndcapRmax']}" deltaphi="360" deltatheta="${attributes['EndcapThetaDeg']}" aunit="deg" lunit="cm"/>

      <union name="SteelVesselUnion1">
       <first ref="SteelTube"/>
       <second ref="EndCap"/>
       <positionref ref="posEndCap1"/>
      </union>
      <union name="SteelVessel">
       <first ref="SteelVesselUnion1"/>
       <second ref="EndCap"/>
       <positionref ref="posEndCap2"/>
       <rotationref ref="rPlus180AboutY"/>
      </union>
      <box name="aSideBeam" lunit="cm" x="256" y="2" z="5"/>
      <box name="aTopBeam" lunit="cm" x="256" y="${4*2.54}" z="2.54"/>
      <box name="aSideCross0" lunit="cm" x="10" y="${0.75*2.54}" z="${43*2.54}"/>
      <box name="aTopCrossBeamA" lunit="cm" x="10" y="${.75*2.54}" z="${49.9*2.54}"/>
      <box name="aTopCrossBeamB" lunit="cm" x="10" y="${.75*2.54}" z="${38*2.54}"/>
      <box name="FrameA" x="11" y="254" z="1040" lunit="cm"/>
      <box name="FrameB" x="11.1" y="${86*2.54}" z="1003.52" lunit="cm"/>
      <box name="VertBar" x="${2.5*2.54}" y="${86*2.54}" z="${2.5*2.54}" lunit="cm"/>
      <box name="CrossBeamA" lunit="cm" x="9" y="301.12" z="7"/>

     <subtraction name="Frame0">
       <first ref="FrameA"/> <second ref="FrameB"/>
       <position name="posFrameSubtraction" x="0" y="0" z="0"/>
     </subtraction>

     <union name="CrossBeam">
        <first ref="CrossBeamA"/> <second ref="CrossBeamA"/>
        <rotation name="CrossBeam93rot" unit="deg" x="93" y="0" z="0"/>
        <position name="posCrossBeamUnion0" unit="cm" x="0" y="0" z="0"/>
      </union>
      <union name="aTopCross0">
        <first ref="aTopCrossBeamA"/> <second ref="aTopCrossBeamA"/>
        <rotationref ref="rPlus60AboutY"/>
        <position name="posaTopCrossUnion" unit="cm" x="0" y="0" z="0"/>
      </union>
      <union name="aTopCross1">
         <first ref="aTopBeam"/> <second ref="aTopCross0"/>
         <rotation name="rMinus60AboutY0" unit="deg" x="0" y="-30" z="0"/>
         <position name="posaCrossBeamUnion0" unit="cm" x="64" y="0" z="${111.64/2}"/>
     </union>
      <union name="aTopCross2">
         <first ref="aTopCross1"/> <second ref="aTopCross0"/>
         <rotation name="rMinus60AboutY1" unit="deg" x="0" y="-30" z="0"/>
         <position name="posaCrossBeamUnion1" unit="cm" x="-64" y="0" z="${111.64/2}"/>
     </union>
     <union name="aTopCross">
        <first ref="aTopCross2"/> <second ref="aTopBeam"/>
        <position name="posaCrossBeamUnion2" unit="cm" x="0" y="0" z="111.64"/>
     </union>

     <union name="aTopCrossOuter0">
       <first ref="aTopCrossBeamB"/> <second ref="aTopCrossBeamB"/>
       <rotation name="rPlus76AboutY" unit="deg" x="0" y="76" z="0" />
       <position name="posaTopCrossOutUnion" unit="cm" x="0" y="0" z="0"/>
     </union>
     <union name="aTopCrossOuter1">
        <first ref="aTopBeam"/> <second ref="aTopCrossOuter0"/>
        <position name="posaCrossOuterUnion0" unit="cm" x="64" y="0" z="36"/>
        <rotation name="rotTopCross1" unit="deg" x="0" y="-38" z="0"/>
    </union>
     <union name="aTopCrossOuter2">
        <first ref="aTopCrossOuter1"/> <second ref="aTopCrossOuter0"/>
        <position name="posaCrossOuterUnion1" unit="cm" x="-64" y="0" z="36"/>
        <rotation name="rotTopCross2" unit="deg" x="0" y="-38" z="0"/>
    </union>
     <union name="aTopCrossOuter">
        <first ref="aTopCrossOuter2"/> <second ref="aTopBeam"/>
        <position name="posaCrossBeamUnion3" unit="cm" x="0" y="0" z="72"/>
    </union>
    <union name="aSideCross1">
      <first ref="aSideCross0"/> <second ref="aSideCross0"/>
      <rotation name="rotSideCross33AboutY" unit="deg" x="0" y="66" z="0"/>
      <position name="posaSideCrossUnion" unit="cm" x="0" y="0" z="0"/>
    </union>
    <union name="aSideCross2">
       <first ref="aTopBeam"/> <second ref="aSideCross1"/>
       <rotation name="rMinus17AboutY" unit="deg" x="0" y="-33" z="0"/>
       <position name="posaSideCrossBeamUnion0" unit="cm" x="64" y="0" z="45.63"/>
   </union>
    <union name="aSideCross3">
       <first ref="aSideCross2"/> <second ref="aSideCross1"/>
       <rotation name="rMinus60AboutY" unit="deg" x="0" y="-33" z="0"/>
       <position name="posaSideCrossBeamUnion1" unit="cm" x="-64" y="0" z="45.63"/>
   </union>
    <union name="aSideCross">
       <first ref="aSideCross2"/> <second ref="aTopBeam"/>
       <position name="posaSideCrossBeamUnion2" unit="cm" x="0" y="0" z="${45.63*2}"/>
   </union>
   %for i in range(1,5):
    <union name="Frame${i}">
       <first ref="Frame${i-1}"/> <second ref="VertBar"/>
       <position name="posVertBarUnion${i}" unit="cm" x="${-11/2.}" y="0" z="${-200.70*3/2. + 200.70*(i-1)}"/>
    </union>
   %endfor

%endif
