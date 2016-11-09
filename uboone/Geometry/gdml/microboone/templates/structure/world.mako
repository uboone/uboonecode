
  <volume name="volVacuumSpace">
    <materialref ref="Vacuum"/>
 	<solidref ref="VacuumSpace"/>
  </volume>

  <volume name="volIBeam22">
	<materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
	<solidref ref="IBeam22"/>
  </volume>
  <volume name="volIBeam35">
	<materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
	<solidref ref="IBeam35"/>
  </volume>

  <volume name="volConcreteWallAboveGrade" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteWallAboveGrade" />
  </volume>
  <volume name="volConcreteDiscRoof" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteDiscRoof" />
  </volume>
  <volume name="volConcreteSteelBeam" >
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
    <solidref ref="ConcreteSteelBeam" />
  </volume>
  <volume name="volSteelGrating" >
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
    <solidref ref="SteelGrating" />
  </volume>
  <volume name="volGround" >
    <materialref ref="Dirt" />
    <solidref ref="Ground" />
  </volume>
  <volume name="volGroundBottom" >
     <materialref ref="Dirt" />
     <solidref ref="GroundBottom" />
   </volume>
  <volume name="volOverburden" >
    <materialref ref="Concrete" />
    <solidref ref="Overburden" />
  </volume>
  <volume name="volPolystyreneEnclosure" >
    <materialref ref="Polystyrene" />
    <solidref ref="PolystyreneEnclosure" />
  </volume>
  <volume name="volConcreteEnclosure" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteEnclosure" />
  </volume>
  <volume name="volPolystyreneEnclosureBottom">
     <materialref ref="Polystyrene" />
     <solidref ref="PolystyreneEnclosureBottom"/>
  </volume>
  <volume name="volConcreteEnclosureBottom" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteEnclosureBottom" />
  </volume>
  <volume name="volRemovableRoof" >
	<materialref ref="Air" />
	<solidref ref="RemovableRoof"/>
 	<physvol>
       <volumeref ref="volIBeam22"/>
       <position name="posIBeam220" unit="cm" x="${(3*12+5)*2.54+4.03*0.5*2.54}" y="0" z="${-(12.5*2.54*0.5)}"/>
     </physvol>
	 <physvol>
       <volumeref ref="volIBeam22"/>
       <position name="posIBeam221" unit="cm" x="${(3*12+5)*2.54+(6*12+10)*2.54+4.03*0.5*2.54}" y="0" z="${-(12.5*2.54*0.5)}"/>
     </physvol>
	 <physvol>
       <volumeref ref="volIBeam22"/>
       <position name="posIBeam222" unit="cm" x="${-((3*12+5)*2.54+4.03*0.5*2.54)}" y="0" z="${-(12.5*2.54*0.5)}"/>
     </physvol>
	 <physvol>
       <volumeref ref="volIBeam22"/>
       <position name="posIBeam223" unit="cm" x="${-((3*12+5)*2.54+(6*12+10)*2.54+4.03*0.5*2.54)}" y="0" z="${-(12.5*2.54*0.5)}"/>
     </physvol>
  %for i in range(4):
		<physvol>
			<volumeref ref="volIBeam35"/>
			<position name="posIBeam35pos${i}" unit="cm" x="0" y="${38*0.5*2.54+i*(6*12+6)*2.54}" z="${(12.5*2.54*0.5)}"/>
			<rotation name="r90AboutZpos${i}" unit="deg" x="0" y="0" z="90"/>
		</physvol>
		<physvol>
			<volumeref ref="volIBeam35"/>
			<position name="posIBeam35neg${i}" unit="cm" x="0" y="${-(38*0.5*2.54+i*(6*12+6)*2.54)}" z="${(12.5*2.54*0.5)}"/>
			<rotation name="r90AboutZneg${i}" unit="deg" x="0" y="0" z="90"/>
		</physvol>

  %endfor

  </volume>
  <volume name="volRoofPlusBars">
	<materialref ref="Air"/>
	<solidref ref="RoofPlusBars"/>
	<physvol>
		<volumeref ref="volSteelGrating"/>
		<position name="posSteelGrating" unit="cm" x="0" y="0" z="50"/>
        <rotationref ref="rPlus90AboutX"/>
	</physvol>
	<physvol>
		<volumeref ref="volConcreteDiscRoof"/>
    	<position name="posConcreteDiscRoof" unit="cm" x="0" y="0" z="0"/>
	</physvol>
	<physvol>
		<volumeref ref="volConcreteSteelBeam"/>
    	<position name="posCBeam1" unit="cm" x="${-6*27*2.54}" y="0" z="0"/>
	</physvol>
	<physvol>
		<volumeref ref="volConcreteSteelBeam"/>
    	<position name="posCBeam2" unit="cm" x="${6*27*2.54}" y="0" z="0"/>
	</physvol>
	<physvol>
		<volumeref ref="volRemovableRoof"/>
		<position name="posRemovableRoof" unit="cm" x="0" y="0" z="0"/>
	</physvol>
	</volume>

  <volume name="volWorld" >
    <materialref ref="Air"/>
    <solidref ref="World"/>
   <physvol>
	  <volumeref ref="volVacuumSpace"/>
 	  <position name="posVacuumSpace" unit="cm" x="0" y="${(53000-1800)/2 + 624.85}" z="0"/>
 	</physvol>
	<physvol>
      <volumeref ref="volConcreteWallAboveGrade"/>
      <position name="posConcreteWallAboveGrade" unit="cm" x="${0.5*256.35}" y="${34.50000555*12*2.54 + 22.86}" z="${0.5*1037}"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
	<physvol>
      <volumeref ref="volRoofPlusBars"/>
      <position name="posRoofPlusBars" unit="cm" x="${0.5*256.35}" y="${34.50000555*12*2.54 +22.86 +127*2.54}" z="${0.5*1037}"/>
	  <rotation name="r39Plus90" unit="deg" x="90" y="0" z="39" />
    </physvol>
    <physvol>
      <volumeref ref="volConcreteEnclosure"/>
      <position name="posConcreteEnclosure" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="${36*2.54/2}" z="${0.5*attributes['TPCWirePlaneLength']}"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
      <volumeref ref="volConcreteEnclosureBottom"/>
      <position name="posConcreteEnclosureBottom" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="${-38*12*2.54/2}" z="${0.5*attributes['TPCWirePlaneLength']}"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
      <volumeref ref="volPolystyreneEnclosure"/>
      <position name="posPolystyreneEnclosure" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="0" z="${0.5*attributes['TPCWirePlaneLength']}"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
       <volumeref ref="volGround"/>
      <position name="posGround" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="0" z="${0.5*attributes['TPCWirePlaneLength']}"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
       <volumeref ref="volGroundBottom"/>
      <position name="posGroundBottom" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="${-41*12*2.54/2 -50*12*2.54/2}" z="${0.5*attributes['TPCWirePlaneLength']}"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
  <!-- <physvol>
      <volumeref ref="volOverburden"/>
      <position name="posOverburden" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="${32.50000555*12*2.54+22.86 +127*2.54+285}" z="${0.5*1037}"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>-->
    <physvol>
      <volumeref ref="volDetEnclosure"/>
      <position name="posDetEnclosure" unit="cm" x="${0.5*attributes['TPCActiveDepth']}" y="0" z="${0.5*attributes['TPCWirePlaneLength']}"/>
    </physvol>
  </volume>
