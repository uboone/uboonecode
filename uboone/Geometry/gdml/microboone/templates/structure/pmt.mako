%if attributes['pmt_on']:
  <volume name="volOpDetSensitive">
   <materialref ref="LAr"/>
   <solidref ref="PMT_AcrylicPlate"/>
  </volume>
  <volume name="vol_PMT_AcrylicPlate">
   <materialref ref="Acrylic"/>
   <solidref ref="PMT_AcrylicPlate"/>
  </volume>
  <volume name="vol_PMT_Stalk">
   <materialref ref="Glass"/>
   <solidref ref="PMT_Stalk"/>
  </volume>
  <volume name="vol_PMT_SteelBase">
   <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <solidref ref="PMT_SteelBase"/>
  </volume>
  <volume name="vol_PMT_Underside">
   <materialref ref="Glass"/>
   <solidref ref="PMT_Underside"/>
  </volume>
  <volume name="volPaddle_PMT">
   <materialref ref="Acrylic"/>
   <solidref ref="Paddle_PMT"/>
  </volume>
  <volume name="vol_PMT_Lens">
   <materialref ref="LAr"/>
   <solidref ref="PMT_Lens"/>
  </volume>
  <volume name="volPMT">
   <materialref ref="LAr"/>
   <solidref ref="PMTVolume"/>
   <physvol>
    <volumeref ref="volOpDetSensitive"/>
    <position name="posOpDetSensitive" unit="cm" x="0" y="0" z="${(5.5 * 2.54) - 0.1}"/>
   </physvol>
   <physvol>
    <volumeref ref="vol_PMT_AcrylicPlate"/>
    <position name="pos_PMT_AcrylicPlate" unit="cm" x="0" y="0" z="${(5.5 * 2.54) - 0.3}"/>
   </physvol>
   <physvol>
    <volumeref ref="vol_PMT_Stalk"/>
    <position name="pos_PMT_Stalk" unit="cm" x="0" y="0" z="${(3.0 * 2.54)-(5.5 * 2.54)}"/>
   </physvol>
   <physvol>
    <volumeref ref="vol_PMT_SteelBase"/>
    <position name="pos_PMT_SteelBase" unit="cm" x="0" y="0" z="${(0.75 * 2.54)-(5.5 * 2.54)}"/>
   </physvol>
   <physvol>
    <volumeref ref="vol_PMT_Lens"/>
    <position name="pos_PMT_Lens" unit="cm" x="0" y="0" z="${(7.0 * 2.54)-(5.5 * 2.54)}"/>
   </physvol>
   <physvol>
    <volumeref ref="vol_PMT_Underside"/>
    <position name="pos_PMT_Underside" unit="cm" x="0" y="0" z="${(7.0 * 2.54)-(5.5 * 2.54)}"/>
   </physvol>
  </volume>
%endif
