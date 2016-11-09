%if attributes['gen_granite']:
<!-- 06-Dec-2011 WGS:
 Cut-and-paste from the GDML description of the Double Chooz Outer Veto module,
 with a change: I assume we don't need multiple definitions of "Air" -->
<volume name="ScintPlastic_3225_log0x1bedfc0">
  <materialref ref="plastic0x1bec080"/>
  <solidref ref="ScintPlastic_32250x1bf1130"/>
</volume>
<volume name="ScintStrip_3225_log0x1beded0">
  <materialref ref="TitaniumDioxide0x1beaf50"/>
  <solidref ref="ScintStrip_32250x1bf10e0"/>
  <physvol name="ScintPlastic_3225_phy0x1bee0b0">
    <volumeref ref="ScintPlastic_3225_log0x1bedfc0"/>
  </physvol>
</volume>
<volume name="TapeSheet_3225_log0x1beeee0">
  <materialref ref="foamtape0x1bec350"/>
  <solidref ref="TapeSheet_32250x1beee50"/>
</volume>
<volume name="SideSpacer_3225_log0x1bef4a0">
  <materialref ref="rubber0x1beccc0"/>
  <solidref ref="SideSpacer_32250x1bef410"/>
</volume>
<volume name="AirBox_log0x1bf0810">
  <materialref ref="Air"/>
  <solidref ref="AirBox0x1bede80"/>
</volume>
<volume name="Module_3225_log0x1bf1050">
  <materialref ref="Aluminum0x1beb1f0"/>
  <solidref ref="Module_32250x1bf1000"/>
  <physvol name="phys_strip_3225_00x1bee160">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_00x1bee160_pos" unit="mm" x="789.075" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_10x1bee1f0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_10x1bee1f0_pos" unit="mm" x="738.975" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_20x1bee280">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_20x1bee280_pos" unit="mm" x="688.875" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_30x1bee340">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_30x1bee340_pos" unit="mm" x="638.775" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_40x1bee3d0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_40x1bee3d0_pos" unit="mm" x="588.675" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_50x1bee4b0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_50x1bee4b0_pos" unit="mm" x="538.575" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_60x1bee540">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_60x1bee540_pos" unit="mm" x="488.475" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_70x1bee5d0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_70x1bee5d0_pos" unit="mm" x="438.375" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_80x1bee660">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_80x1bee660_pos" unit="mm" x="388.275" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_90x1bee420">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_90x1bee420_pos" unit="mm" x="338.175" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_100x1b49940">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_100x1b49940_pos" unit="mm" x="288.075" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_110x1b499d0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_110x1b499d0_pos" unit="mm" x="237.975" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_120x1b49a60">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_120x1b49a60_pos" unit="mm" x="187.875" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_130x1b49af0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_130x1b49af0_pos" unit="mm" x="137.775" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_140x1b49b80">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_140x1b49b80_pos" unit="mm" x="87.675" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_150x1b49c10">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_150x1b49c10_pos" unit="mm" x="37.575" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_160x1b49ca0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_160x1b49ca0_pos" unit="mm" x="-12.525" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_170x1b49870">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_170x1b49870_pos" unit="mm" x="-62.625" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_180x1b49e40">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_180x1b49e40_pos" unit="mm" x="-112.725" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_190x1b49ed0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_190x1b49ed0_pos" unit="mm" x="-162.825" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_200x1b49f60">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_200x1b49f60_pos" unit="mm" x="-212.925" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_210x1b49ff0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_210x1b49ff0_pos" unit="mm" x="-263.025" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_220x1b4a080">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_220x1b4a080_pos" unit="mm" x="-313.125" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_230x1b4a110">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_230x1b4a110_pos" unit="mm" x="-363.225" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_240x1b4a1a0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_240x1b4a1a0_pos" unit="mm" x="-413.325" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_250x1b4a230">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_250x1b4a230_pos" unit="mm" x="-463.425" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_260x1b4a2c0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_260x1b4a2c0_pos" unit="mm" x="-513.525" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_270x1b4a350">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_270x1b4a350_pos" unit="mm" x="-563.625" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_280x1b4a3e0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_280x1b4a3e0_pos" unit="mm" x="-613.725" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_290x1b4a470">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_290x1b4a470_pos" unit="mm" x="-663.825" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_300x1b4a500">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_300x1b4a500_pos" unit="mm" x="-713.925" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_310x1b4a590">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_310x1b4a590_pos" unit="mm" x="-764.025" y="-6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_320x1b4a620">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_320x1b4a620_pos" unit="mm" x="764.025" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_330x1b49d30">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_330x1b49d30_pos" unit="mm" x="713.925" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_340x1b4a880">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_340x1b4a880_pos" unit="mm" x="663.825" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_350x1b4a8d0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_350x1b4a8d0_pos" unit="mm" x="613.725" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_360x1b4a960">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_360x1b4a960_pos" unit="mm" x="563.625" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_370x1b4a9f0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_370x1b4a9f0_pos" unit="mm" x="513.525" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_380x1b4aa80">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_380x1b4aa80_pos" unit="mm" x="463.425" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_390x1b4ab10">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_390x1b4ab10_pos" unit="mm" x="413.325" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_400x1b4aba0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_400x1b4aba0_pos" unit="mm" x="363.225" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_410x1b4ac30">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_410x1b4ac30_pos" unit="mm" x="313.125" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_420x1b4acc0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_420x1b4acc0_pos" unit="mm" x="263.025" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_430x1b4ad50">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_430x1b4ad50_pos" unit="mm" x="212.925" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_440x1b4ade0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_440x1b4ade0_pos" unit="mm" x="162.825" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_450x1b4ae70">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_450x1b4ae70_pos" unit="mm" x="112.725" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_460x1b4af00">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_460x1b4af00_pos" unit="mm" x="62.625" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_470x1b4af90">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_470x1b4af90_pos" unit="mm" x="12.525" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_480x1b4b020">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_480x1b4b020_pos" unit="mm" x="-37.575" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_490x1b4b0b0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_490x1b4b0b0_pos" unit="mm" x="-87.675" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_500x1bee6b0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_500x1bee6b0_pos" unit="mm" x="-137.775" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_510x1bee740">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_510x1bee740_pos" unit="mm" x="-187.875" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_520x1bee7d0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_520x1bee7d0_pos" unit="mm" x="-237.975" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_530x1bee860">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_530x1bee860_pos" unit="mm" x="-288.075" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_540x1bee8f0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_540x1bee8f0_pos" unit="mm" x="-338.175" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_550x1bee980">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_550x1bee980_pos" unit="mm" x="-388.275" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_560x1beea10">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_560x1beea10_pos" unit="mm" x="-438.375" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_570x1beeaa0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_570x1beeaa0_pos" unit="mm" x="-488.475" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_580x1beeb30">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_580x1beeb30_pos" unit="mm" x="-538.575" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_590x1beebc0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_590x1beebc0_pos" unit="mm" x="-588.675" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_600x1beec50">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_600x1beec50_pos" unit="mm" x="-638.775" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_610x1beece0">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_610x1beece0_pos" unit="mm" x="-688.875" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_620x1beed70">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_620x1beed70_pos" unit="mm" x="-738.975" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_strip_3225_630x1beee00">
    <volumeref ref="ScintStrip_3225_log0x1beded0"/>
    <position name="phys_strip_3225_630x1beee00_pos" unit="mm" x="-789.075" y="6.2032" z="190"/>
  </physvol>
  <physvol name="phys_tapesheet_3225_00x1beefb0">
    <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
    <position name="phys_tapesheet_3225_00x1beefb0_pos" unit="mm" x="0" y="11.7032" z="190"/>
  </physvol>
  <physvol name="phys_tapesheet_3225_10x1b4a6b0">
    <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
    <position name="phys_tapesheet_3225_10x1b4a6b0_pos" unit="mm" x="0" y="0.7032" z="190"/>
  </physvol>
  <physvol name="phys_tapesheet_3225_20x1b4a740">
    <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
    <position name="phys_tapesheet_3225_20x1b4a740_pos" unit="mm" x="0" y="-0.7032" z="190"/>
  </physvol>
  <physvol name="phys_tapesheet_3225_30x1b4a7d0">
    <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
    <position name="phys_tapesheet_3225_30x1b4a7d0_pos" unit="mm" x="0" y="-11.7032" z="190"/>
  </physvol>
  <physvol name="phys_spacer3225_10x1bef570">
    <volumeref ref="SideSpacer_3225_log0x1bef4a0"/>
    <position name="phys_spacer3225_10x1bef570_pos" unit="mm" x="819.125" y="0" z="0"/>
  </physvol>
  <physvol name="phys_spacer3225_20x1bef600">
    <volumeref ref="SideSpacer_3225_log0x1bef4a0"/>
    <position name="phys_spacer3225_20x1bef600_pos" unit="mm" x="-819.125" y="0" z="0"/>
  </physvol>
  <physvol name="phys_airbox0x1bef650">
    <volumeref ref="AirBox_log0x1bf0810"/>
    <position name="phys_airbox0x1bef650_pos" unit="mm" x="0" y="0" z="-1612.6"/>
  </physvol>
</volume>

<!-- 06-Dec-11 WGS:
  To rotate the scintillator modules properly, we need two separate
  rotations. Take this opportunity to define the scintillator modules with
  four distinct names. We'll need to identify the volumes separately when
  we write hit-processing code for Geant4. -->
<volume name="volDoubleChoozScintModule_0">
  <materialref ref="Aluminum0x1beb1f0"/>
  <solidref ref="Module_3225_rot"/>
  <physvol>
    <volumeref ref="Module_3225_log0x1bf1050" />
        <positionref ref="posCenter" />
        <rotationref ref="rPlus90AboutZ" />
  </physvol>
</volume>
<volume name="volDoubleChoozScintModule_1">
  <materialref ref="Aluminum0x1beb1f0"/>
  <solidref ref="Module_3225_rot"/>
  <physvol>
    <volumeref ref="Module_3225_log0x1bf1050" />
        <positionref ref="posCenter" />
        <rotationref ref="rPlus90AboutZ" />
  </physvol>
</volume>
<volume name="volDoubleChoozScintModule_2">
  <materialref ref="Aluminum0x1beb1f0"/>
  <solidref ref="Module_3225_rot"/>
  <physvol>
    <volumeref ref="Module_3225_log0x1bf1050" />
        <positionref ref="posCenter" />
        <rotationref ref="rPlus90AboutZ" />
  </physvol>
</volume>
<volume name="volDoubleChoozScintModule_3">
  <materialref ref="Aluminum0x1beb1f0"/>
  <solidref ref="Module_3225_rot"/>
  <physvol>
    <volumeref ref="Module_3225_log0x1bf1050" />
        <positionref ref="posCenter" />
        <rotationref ref="rPlus90AboutZ" />
  </physvol>
</volume>
<volume name="volGraniteBlock">
 <materialref ref="Granite" />
 <solidref ref="GraniteBlock" />
</volume>
%endif
