%if attributes['gen_extras']:

  <volume name="volInsulation">
     <materialref ref="PU_foam_light"/>
     <solidref ref="Insulation"/>
  </volume>
  <volume name="volPlatformDisc">
  <materialref ref="G10"/>
    <solidref ref="PlatformDisc"/>
  </volume>
  <volume name="volpolySmallSquare">
     <materialref ref="G10"/>
     <solidref ref="polySmallSquare"/>
  </volume>
  <volume name="volSquare">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="Square"/>
  </volume>
  <volume name="volIBeam53">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="IBeam53"/>
  </volume>
  <volume name="volIBeam16">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="IBeam16"/>
  </volume>
  <volume name="volIBeam40">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="IBeam40"/>
  </volume>
  <volume name="volIBeam17">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="IBeam17"/>
  </volume>
  <volume name="volIBeam45">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="IBeam45"/>
  </volume>

    <volume name="volPlatformSpace">
        <materialref ref="Air"/>
        <solidref ref="PlatformSpace"/>
  <physvol>
        <volumeref ref="volPlatformDisc"/>
        <position name="posPlatformDisc" unit="cm" x="0" y="${(4*12+8+ 11.94*2+2)*2.54/2}" z="0"/>
      </physvol>
  %for i in range(5):
      <physvol>
        <volumeref ref="volIBeam40"/>
        <position name="posIbeam40pos${i}" unit="cm" x="0" y="${(4*12+8-11.94)*2.54/2-.5}" z="${(24.5+49*i)*2.54}"/>
        <rotation name="rPlus90Z90X40pos${i}" unit="deg" x="90" y="0" z="90"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam40"/>
        <position name="posIbeam40neg${i}" unit="cm" x="0" y="${(4*12+8-11.94)*2.54/2-.5}" z="${(-24.5-49*i)*2.54}"/>
        <rotation name="rPlus90Z90X40neg${i}" unit="deg" x="90" y="0" z="90"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam17"/>
        <position name="posIbeam17pos${i}" unit="cm" x="${(-9*12-2+18)*2.54+.5}" y="${-(4*12+8-10.11)*2.54/2+.5-10.11*2.54}" z="${(24.5+49*i)*2.54}"/>
        <rotation name="rPlus90Z90Xpos17${i}" unit="deg" x="90" y="0" z="90"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam17"/>
        <position name="posIbeam17neg${i}" unit="cm" x="${(-9*12-2+18)*2.54+.5}" y="${-(4*12+8-10.11)*2.54/2+.5-10.11*2.54}" z="${(-24.5-49*i)*2.54}"/>
        <rotation name="rPlus90Z90Xneg17${i}" unit="deg" x="90" y="0" z="90"/>
      </physvol>
  %endfor
  %for j in range(3):
      <physvol>
        <volumeref ref="volSquare"/>
        <position name="posSquareposZ${j}" unit="cm" x="${(-9*12-2+4)*2.54+.5}" y="${-(4*12+8-42)*2.54/2+.5}" z="${(24.5+2*49*j)*2.54}"/>
      </physvol>
      <physvol>
        <volumeref ref="volSquare"/>
        <position name="posSquarenegZ${j}" unit="cm" x="${(-9*12-2+4)*2.54+.5}" y="${-(4*12+8-42)*2.54/2+.5}" z="${-(24.5+2*49*j)*2.54}"/>
      </physvol>
      <physvol>
        <volumeref ref="volSquare"/>
        <position name="posSquareposZShift${j}" unit="cm" x="${(-9*12-2+8+36)*2.54+1.25}" y="${-(4*12+8-42)*2.54/2+.5}" z="${(24.5+2*49*j)*2.54}"/>
      </physvol>
      <physvol>
        <volumeref ref="volSquare"/>
        <position name="posSquarenegZShift${j}" unit="cm" x="${(-9*12-2+8+36)*2.54+1.25}" y="${-(4*12+8-42)*2.54/2+.5}" z="${-(24.5+2*49*j)*2.54}"/>
      </physvol>
  %endfor
  %for j in range(2):
      <physvol>
        <volumeref ref="volIBeam53"/>
        <position name="posIBeam53pos${j}" unit="cm" x="${(-1)**j*((9*12+2-9.995/2)*2.54-.5)}" y="${(4*12+8+12.06)*2.54/2}" z="${(23*12/2-5)*2.54}"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam53"/>
        <position name="posIBeam53neg${j}" unit="cm" x="${(-1)**j*((9*12+2-9.995/2)*2.54-.5)}" y="${(4*12+8+12.06)*2.54/2}" z="${-(23*12/2-5)*2.54}"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam16"/>
        <position name="posIBeam16pos${j}" unit="cm" x="${(-1)**j*((6*12)/2*2.54-.5)}" y="${(4*12+8+12.06)*2.54/2}" z="${128.5*2.54}"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam16"/>
        <position name="posIBeam16neg${j}" unit="cm" x="${(-1)**j*((6*12)/2*2.54-.5)}" y="${(4*12+8+12.06)*2.54/2}" z="${-128.5*2.54}"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam45"/>
        <position name="posIBeam45pos${j}" unit="cm" x="${-(9*12+2-8-8.02/2-.5)*2.54 +j*(36-1.5*8.02)*2.54 }" y="${-(4*12+8-2*10.11)*2.54/2}" z="${128.5*2.54}"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volIBeam45"/>
        <position name="posIBeam45neg${j}" unit="cm" x="${-(9*12+2-8-8.02/2-.5)*2.54+j*(36-1.5*8.02)*2.54}" y="${-(4*12+8-2*10.11)*2.54/2}" z="${-128.5*2.54}"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
  %endfor
      <physvol>
        <volumeref ref="volpolySmallSquare"/>
        <position name="pospolySmallSquare" unit="cm" x="${-(9*12+2-18-10.11/2-.5)*2.54}" y="${-(4*12+8-3*10.11-3)*2.54/2}" z="${0}"/>
      </physvol>

    </volume>
    <volume name="volColumn">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="Column"/>
    </volume>
    <volume name="volTankCap1">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="TankCap1"/>
    </volume>
    <volume name="volTank1">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="Tank1"/>
    </volume>
    <volume name="volStandSubtraction">
      <materialref ref="PU_foam_dense"/>
      <solidref ref="standSubtraction"/>
    </volume>

    <volume name="volTankBox1">
      <materialref ref="Air"/>
      <solidref ref="TankBox1"/>
      <physvol>
        <volumeref ref="volTank1"/>
        <position name="posvolTank1" unit="cm" x="0" y="0" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volTankCap1"/>
        <position name="posvolTankCap1_1" unit="cm" x="0" y="${90.15+.001}" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volTankCap1"/>
        <position name="posvolTankCap1_2" unit="cm" x="0" y="${-90.15-.001}" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
    </volume>

    <volume name="volStandConcrete">
      <materialref ref="Concrete"/>
      <solidref ref="standConcrete"/>
    </volume>

    <volume name="volRackX">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="rackX"/>
    </volume>
    <volume name="volRackY">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="rackY"/>
    </volume>
    <volume name="volRackZ">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="rackZ"/>
    </volume>
  <volume name="volRackInnards">
    <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
    <solidref ref="rackInnards"/>
  </volume>

    <volume name="volRack">
      <materialref ref="Air"/>
      <solidref ref="rackBox"/>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX1" unit="cm" x="0" y="${(attributes['RackY']-attributes['RackThickness'])/2}" z="${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX2" unit="cm" x="0" y="-${(attributes['RackY']-attributes['RackThickness'])/2}" z="${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX3" unit="cm" x="0" y="${(attributes['RackY']-attributes['RackThickness'])/2}" z="-${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX4" unit="cm" x="0" y="-${(attributes['RackY']-attributes['RackThickness'])/2}" z="-${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY1" unit="cm" y="0" x="${(attributes['RackX']-attributes['RackThickness'])/2}" z="${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY2" unit="cm" y="0" x="-${(attributes['RackX']-attributes['RackThickness'])/2}" z="${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY3" unit="cm" y="0" x="${(attributes['RackX']-attributes['RackThickness'])/2}" z="-${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY4" unit="cm" y="0" x="-${(attributes['RackX']-attributes['RackThickness'])/2}" z="-${(attributes['RackZ']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ1" unit="cm" z="0" y="${(attributes['RackY']-attributes['RackThickness'])/2}" x="${(attributes['RackX']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ2" unit="cm" z="0" y="${(attributes['RackY']-attributes['RackThickness'])/2}" x="-${(attributes['RackX']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ3" unit="cm" z="0" y="-${(attributes['RackY']-attributes['RackThickness'])/2}" x="${(attributes['RackX']-attributes['RackThickness'])/2}"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ4" unit="cm" z="0" y="-${(attributes['RackY']-attributes['RackThickness'])/2}" x="-${(attributes['RackX']-attributes['RackThickness'])/2}"/>
      </physvol>
    <physvol>
    <volumeref ref="volRackInnards"/>
    <position name="posRackInnards" unit="cm" x="0" y="0" z="0"/>
    </physvol>
    </volume>

  <volume name="volFloorTank1">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="floorTank1"/>
  </volume>

  <volume name="volFloorTankCap1">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="floorTankCap1"/>
  </volume>

  <volume name="volFloorTank2">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="floorTank2"/>
  </volume>

  <volume name="volFloorTankCap2">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="floorTankCap2"/>
  </volume>

  <volume name="volFloorTankBox1">
      <materialref ref="Air"/>
      <solidref ref="floorTankBox1"/>
      <physvol>
        <volumeref ref="volFloorTank1"/>
        <position name="posFloorTank1" unit="cm" x="0" y="0" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volFloorTankCap1"/>
        <position name="posFloorTankCap1_1" unit="cm" x="0" y="89.9001" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volFloorTankCap1"/>
        <position name="posFloorTankCap1_2" unit="cm" x="0" y="-89.9001" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
    </volume>

    <volume name="volFloorTankBox2">
      <materialref ref="Air"/>
      <solidref ref="floorTankBox2"/>
      <physvol>
        <volumeref ref="volFloorTank2"/>
        <position name="posFloorTank2" unit="cm" x="0" y="6.033" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volFloorTankCap2"/>
        <position name="posFloorTankCap2" unit="cm" x="0" y="-26.98" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
    </volume>

     <volume name="volWalkway">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="Walkway"/>
     </volume>


     <volume name="volPump">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="Pump"/>
     </volume>


%endif
