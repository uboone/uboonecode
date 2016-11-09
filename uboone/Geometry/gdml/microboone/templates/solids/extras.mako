%if attributes['gen_extras']:
  <tube name="Insulation"
    rmin="${attributes['CryostatOuterRadius']+0.1}"
    rmax="${attributes['CryostatOuterRadius']+0.1+attributes['inch']*attributes['InsulationThickness']}"
    z="${attributes['CryostatLength']}"
    deltaphi="360"
    aunit="deg"
    lunit="cm"/>
   <box name="PlatformOLDbig" lunit="cm" x="${558+1}" y="${(4*12+8+2*11+12.06)*2.54+1}" z="${1500+1}" />
   <box name="PlatformOLDsmall" lunit="cm" x="558" y="${(4*12+8+2*11+12.06)*2.54}" z="1500" />
   <tube name="PlatformNEW" lunit="cm" rmax="${292*2.54-0.5}" z="${(4*12+7+2*11+12.06)*2.54+1.5}" aunit="deg" deltaphi="360" />

   <box name="PlatformRectbig" lunit="cm" x="${550+1}" y="${2.54+1}" z="${1500+1}" />
   <box name="PlatformRectsmall" lunit="cm" x="550" y="2.54" z="1500" />
   <tube name="PlatformCircle" lunit="cm" rmax="${292*2.54-0.5}" z="${2.54+1.5}" aunit="deg" deltaphi="360" />

    <subtraction name="PlatformSub0">
      <first ref="PlatformOLDbig"/> <second ref="PlatformNEW"/>
  	<position name="posPlatformSub0" unit="cm" x="0" y="0" z="0"/>
  	<rotation name="rotPlatformSub0" unit="deg" x="90" y="0" z="0"/>
    </subtraction>

    <subtraction name="PlatformSpace0">
  	<first ref="PlatformOLDsmall"/> <second ref="PlatformSub0"/>
  	<position name="posPlatformSub1" unit="cm" x="0" y="0" z="0"/>
    </subtraction>

    <subtraction name="PlatformDisc0">
      <first ref="PlatformRectbig"/> <second ref="PlatformCircle"/>
  	<position name="posPlatformDiscSub" unit="cm" x="0" y="0" z="0"/>
  	<rotation name="rotPlatformDisc" unit="deg" x="90" y="0" z="0"/>
    </subtraction>

    <subtraction name="PlatformDisc">
      <first ref="PlatformRectsmall"/> <second ref="PlatformDisc0"/>
  	<position name="posPlatformDiscSub1" unit="cm" x="0" y="0" z="0"/>
    </subtraction>


    <subtraction name="PlatformSpace">
  	<first ref="PlatformSpace0"/> <second ref="PlatformOLDbig"/>
  	<position name="posPlatformSub2" unit="cm" x="${(45+9)*2.54}" y="${-30*2.54}" z="0"/>
    </subtraction>

   <box name="IBeam53_0" lunit="cm" x="${9.995*2.54 -0.1}" y="${(22*12+6)*2.54-0.1}" z="${12.06*2.54-0.1}"/>
   <box name="IBeam53Segment" lunit="cm" x="${(9.995-0.345)*2.54/2}" y="${(22*12+6)*2.54}" z="${(12.06-2*0.575)*2.54}"/>
   <box name="IBeam16_0" lunit="cm" x="${3.99*2.54-0.1}" y="${(17*12+4)*2.54-0.1}" z="${12.5*2.54-0.1}"/>
   <box name="IBeam16Segment" lunit="cm" x="${(3.99-0.22)*2.54/2}" y="${(17*12+4)*2.54}" z="${(11.99-2*0.265)*2.54}"/>
   <box name="IBeam40_0" lunit="cm" x="${8.005*2.54 -0.1}" y="${(18*12+3)*2.54-0.1}" z="${11.94*2.54-0.1}"/>
   <box name="IBeam40Segment" lunit="cm" x="${(8.005-0.295)*2.54/2}" y="${(18*12+3)*2.54}" z="${(11.94-2*0.515)*2.54}"/>
   <box name="IBeam45_0" lunit="cm" x="${8.02*2.54-0.1}" y="${(17*12 +4)*2.54-0.1}" z="${10.1*2.54-0.1}"/>
   <box name="IBeam45Segment" lunit="cm" x="${(8.02-0.35)*2.54/2}" y="${(17*12+4)*2.54}" z="${(10.1-2*0.62)*2.54}"/>
   <box name="IBeam17_0" lunit="cm" x="${4.01*2.54-0.1}" y="${36*2.54-0.1}" z="${10.11*2.54-0.1}"/>
   <box name="IBeam17Segment" lunit="cm" x="${(4.01-0.24)*2.54/2}" y="${36*2.54}" z="${(10.11-2*0.33)*2.54}"/>
   <box name="BigSquare" lunit="cm" x="${8*2.54}" y="${(4*12+6-12)*2.54}" z="${8*2.54}"/>
   <box name="SmallSquare" lunit="cm" x="${7*2.54}" y="${(4*12+6-10)*2.54}" z="${7*2.54}"/>
   <box name="polySquare" lunit="cm" x="${(18*12+4)*2.54-1}" y="${2.54}" z="${44*12*2.54}"/>
   <box name="polySmallSquare" lunit="cm" x="${(30)*2.54}" y="${2.54}" z="${465*2.54}"/>

    <subtraction name="Square">
      <first ref="BigSquare"/> <second ref="SmallSquare"/>
  	<position name="posSquareRemoval0" unit="cm" x="0" y="0" z="0"/>
    </subtraction>
    <subtraction name="IBeam53_1">
      <first ref="IBeam53_0"/> <second ref="IBeam53Segment"/>
  	<position name="posBeamRemovalPlat2" unit="cm" x="${(9.995+.345)*0.25*2.54}" y="0" z="0"/>
    </subtraction>
    <subtraction name="IBeam53">
      <first ref="IBeam53_1"/> <second ref="IBeam53Segment"/>
  	<position name="posBeamRemovalPlat3" unit="cm" x="${-(9.995 +.345)*0.25*2.54}" y="0" z="0"/>
    </subtraction>

    <subtraction name="IBeam16_1">
      <first ref="IBeam16_0"/> <second ref="IBeam16Segment"/>
  	<position name="posBeamRemovalPlat4" unit="cm" x="${(3.99+.22)*0.25*2.54}" y="0" z="0"/>
    </subtraction>
    <subtraction name="IBeam16">
      <first ref="IBeam16_1"/> <second ref="IBeam16Segment"/>
  	<position name="posBeamRemovalPlat5" unit="cm" x="${-(3.99+.22)*0.25*2.54}" y="0" z="0"/>
    </subtraction>

    <subtraction name="IBeam40_1">
      <first ref="IBeam40_0"/> <second ref="IBeam40Segment"/>
  	<position name="posBeamRemovalPlat6" unit="cm" x="${(8.005+.295)*0.25*2.54}" y="0" z="0"/>
    </subtraction>
    <subtraction name="IBeam40">
      <first ref="IBeam40_1"/> <second ref="IBeam40Segment"/>
  	<position name="posBeamRemovalPlat7" unit="cm" x="${-(8.005+0.295)*0.25*2.54}" y="0" z="0"/>
    </subtraction>

    <subtraction name="IBeam17_1">
      <first ref="IBeam17_0"/> <second ref="IBeam17Segment"/>
  	<position name="posBeamRemovalPlat8" unit="cm" x="${(4.01+.24)*0.25*2.54}" y="0" z="0"/>
    </subtraction>
    <subtraction name="IBeam17">
      <first ref="IBeam17_1"/> <second ref="IBeam17Segment"/>
  	<position name="posBeamRemovalPlat9" unit="cm" x="${-(4.01+.24)*0.25*2.54}" y="0" z="0"/>
    </subtraction>

    <subtraction name="IBeam45_1">
      <first ref="IBeam45_0"/> <second ref="IBeam45Segment"/>
  	<position name="posBeamRemovalPlat10" unit="cm" x="${(8.02+.35)*0.25*2.54}" y="0" z="0"/>
    </subtraction>
    <subtraction name="IBeam45">
      <first ref="IBeam45_1"/> <second ref="IBeam45Segment"/>
  	<position name="posBeamRemovalPlat11" unit="cm" x="${-(8.02+.35)*0.25*2.54}" y="0" z="0"/>
    </subtraction>




  <box name="Column"
    x="16.79"
    y="750"
    z="16.79"
    lunit="cm"/>
  <tube name="Tank1"
    rmin="30.35"
    rmax="38"
    z="140.3"
    deltaphi="360"
    aunit="deg"
    lunit="cm"/>
  <tube name="TankCap1"
    rmax="38"
    z="40"
    deltaphi="360"
    aunit="deg"
    lunit="cm"/>
  <box name="TankBox1"
    x="76.001"
    y="221"
    z="76.001"
    lunit="cm"/>

  <box name="standBox"
    x="507"
    y="196.27"
    z="91.44"
    lunit="cm"/>
  <tube name="standSubTube"
    rmax="233.75"
    z="125"
    deltaphi="360"
    aunit="deg"
    lunit="cm"/>
  <subtraction name="standSubtraction">
    <first ref="standBox"/>
    <second ref="standSubTube"/>
    <position name="posStandSubtract" unit="cm" x="0" y="235.92" z="0"/>
  </subtraction>
  <box name="standConcrete"
    x="548.64"
    y="182.89"
    z="91.44"
    lunit="cm"/>

  <box name="rackBox"
    x="${attributes['RackX']}"
    y="${attributes['RackY']}"
    z="${attributes['RackZ']}"
    lunit="cm"/>
  <box name="rackX"
    x="${attributes['RackX']}"
    y="${attributes['RackThickness']}"
    z="${attributes['RackThickness']}"
    lunit="cm"/>
  <box name="rackY"
    x="${attributes['RackThickness']}"
    y="${attributes['RackY']-2*attributes['RackThickness']-0.001}"
    z="${attributes['RackThickness']}"
    lunit="cm"/>
  <box name="rackZ"
    x="${attributes['RackThickness']}"
    y="${attributes['RackThickness']}"
    z="${attributes['RackZ']-2*attributes['RackThickness']-0.001}"
    lunit="cm"/>

  <box name="rackInnards"
    x="33.0"
    y="84.5"
    z="20.5"
    lunit="cm"/>

  <box name="floorTankBox1"
     x="59.7"
     y="200"
     z="59.7"
     lunit="cm"/>
  <tube name="floorTank1"
     rmin="${45.4/2}"
     rmax="${59.6/2}"
     z="159.8"
     deltaphi="360"
     aunit="deg"
     lunit="cm"/>
  <tube name="floorTankCap1"
     rmax="${59.6/2}"
     z="20.0"
     deltaphi="360"
     aunit="deg"
     lunit="cm"/>
  <box name="floorTankBox2"
     x="52.2"
     y="67"
     z="52.2"
     lunit="cm"/>
  <tube name="floorTank2"
     rmin="${40.0/2}"
     rmax="${52.1/2}"
     z="53.0"
     deltaphi="360"
     aunit="deg"
     lunit="cm"/>
  <tube name="floorTankCap2"
     rmax="${52.1/2}"
     z="${52.1/4}"
     deltaphi="360"
     aunit="deg"
     lunit="cm"/>

  <box name="Walkway"
     x="91"
     y="24"
     z="1200"
     lunit="cm"/>

  <box name="Pump"
     x="100"
     y="196"
     z="91.7"
     lunit="cm"/>

%endif
