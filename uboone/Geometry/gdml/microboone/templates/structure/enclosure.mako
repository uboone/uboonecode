<volume name="volDetEnclosure">
 <materialref ref="Air"/>
 <solidref ref="DetEnclosure"/>
 <physvol>
  <volumeref ref="volCryostat"/>
  <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
 </physvol>
 %if attributes['gen_extras']:
   <physvol>
      <volumeref ref="volInsulation"/>
      <position name="posInsulation" unit="cm" x="0" y="0" z="0"/>
    </physvol>
    <physvol>
      <volumeref ref="volPlatformSpace"/>
      <position name="posPlatform" unit="cm" x="0" y="292.74" z="0"/>
    </physvol>
    <physvol>
      <volumeref ref="volTankBox1"/>
      <position name="posTank1_1" unit="cm" x="50" y="${419+98.1}" z="-600"/>
    </physvol>
    <physvol>
     <volumeref ref="volTankBox1"/>
     <position name="posTank1_2" unit="cm" x="-50" y="${419+98.1}" z="-600"/>
    </physvol>
    <physvol>
     <volumeref ref="volStandSubtraction"/>
     <position name="posStand1" unit="cm" x="0" y="-236" z="335.28"/>
    </physvol>
    <physvol>
       <volumeref ref="volStandSubtraction"/>
     <position name="posStand2" unit="cm" x="0" y="-236" z="-335.28"/>
    </physvol>
    <physvol>
     <volumeref ref="volStandConcrete"/>
     <position name="posConcreteStand1" unit="cm" x="0" y="-424" z="335.28"/>
    </physvol>
    <physvol>
       <volumeref ref="volStandConcrete"/>
     <position name="posConcreteStand2" unit="cm" x="0" y="-424" z="-335.28"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack1" unit="cm" x="-75" y="${408.91+104.1}" z="${0+100+20}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack2" unit="cm" x="-75" y="${408.91+104.1}" z="${93.85+100+20}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack3" unit="cm" x="-75" y="${408.91+104.1}" z="${187.7+100+90}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack4" unit="cm" x="-75" y="${408.91+104.1}" z="${275.15+100+70}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack5" unit="cm" x="-75" y="${408.91+104.1}" z="${-142.15+100}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack6" unit="cm" x="-75" y="${408.91+104.1}" z="${-243.8+100}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack7" unit="cm" x="-75" y="${408.91+104.1}" z="${-326.9+100}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack8" unit="cm" x="127.9" y="${408.91+104.1}" z="${-326.9+122}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack10" unit="cm" x="127.9" y="${408.91+104.1}" z="${-326.9+150+40}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack11" unit="cm" x="100" y="${408.91+104.1}" z="${-387.7-40}"/>
      <rotationref ref="rPlus90AboutY"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack12" unit="cm" x="105" y="${408.91+104.1}" z="${387.7+150+60.8}"/>
      <rotationref ref="rPlus90AboutY"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack13" unit="cm" x="40" y="${408.91+104.1}" z="${387.7+150+60.8}"/>
      <rotationref ref="rPlus90AboutY"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack14" unit="cm" x="-75" y="${408.91+104.1}" z="${270+40}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack15" unit="cm" x="145" y="${408.91+104.1}" z="${-142.15+100}"/>
    </physvol>
    <physvol>
      <volumeref ref="volRack"/>
      <position name="posRack16" unit="cm" x="145" y="${408.91+104.1}" z="170"/>
    </physvol>
    <physvol>
       <volumeref ref="volFloorTankBox1"/>
       <position name="posfloortankbox1" unit="cm" x="-450" y="${100+.001-530}" z="335.28"/>
    </physvol>
    <physvol>
       <volumeref ref="volFloorTankBox1"/>
       <position name="posfloortankbox1_2" unit="cm" x="100" y="${100+.001-530}" z="-445"/>
    </physvol>
    <physvol>
       <volumeref ref="volFloorTankBox1"/>
       <position name="posfloortankbox1_3" unit="cm" x="100" y="${100+.001-530}" z="-521.5"/>
    </physvol>
    <physvol>
       <volumeref ref="volFloorTankBox1"/>
       <position name="posfloortankbox1_4" unit="cm" x="360" y="${100+.001-530}" z="-335.28"/>
    </physvol>
    <physvol>
       <volumeref ref="volFloorTankBox1"/>
       <position name="posfloortankbox1_5" unit="cm" x="${360+60}" y="${100+.001-530}" z="${-335.28-60}"/>
    </physvol>
    <physvol>
       <volumeref ref="volFloorTankBox2"/>
       <position name="posfloortankbox2_1" unit="cm" x="-20" y="${100+.001-530}" z="-500"/>
    </physvol>
    <physvol>
       <volumeref ref="volFloorTankBox2"/>
       <position name="posfloortankbox2_2" unit="cm" x="-136" y="${40+80-530}" z="-500"/>
    </physvol>
    <physvol>
       <volumeref ref="volPump"/>
       <position name="posPump" unit="cm" x="0" y="-424" z="200"/>
    </physvol>
 %endif
 %if attributes['gen_granite']:
 <physvol>
   <volumeref ref="volGraniteBlock"/>
   <position name="posGraniteBlock" unit="cm" x="${-(193+(160/2)+10)}" y="0" z="233"/>
 </physvol>
 <physvol>
   <volumeref ref="volDoubleChoozScintModule_0"/>
   <position name="posDoubleChoozScintModule_0" unit="cm" x="${-(193+(5*0.5))}" y="0" z="233" />
   <rotation name="rotDoubleChoozScintModule_0" unit="deg" x="45" y="0" z="0" />
 </physvol>
 <physvol>
   <volumeref ref="volDoubleChoozScintModule_1"/>
   <position name="posDoubleChoozScintModule_1" unit="cm" x="${-(193+(5*1.5))}" y="0" z="233" />
   <rotation name="rotDoubleChoozScintModule_1" unit="deg" x="135" y="0" z="0" />
 </physvol>
 <physvol>
   <volumeref ref="volDoubleChoozScintModule_2"/>
   <position name="posDoubleChoozScintModule_2" unit="cm" x="${-(193+160+(5*2.5))}" y="0" z="233" />
   <rotation name="rotDoubleChoozScintModule_2" unit="deg" x="45" y="0" z="0" />
 </physvol>
 <physvol>
   <volumeref ref="volDoubleChoozScintModule_3"/>
   <position name="posDoubleChoozScintModule_3" unit="cm" x="${-(193+160+(5*3.5))}" y="0" z="233" />
   <rotation name="rotDoubleChoozScintModule_3" unit="deg" x="135" y="0" z="0" />
 </physvol>
 %endif
 
 %if attributes['gen_vetowall']:
   <%include file="vetowall_placement.mako"/>
 %endif
 %if attributes['gen_crt_a'] and not attributes['gen_crt_b']:
   <%include file="crt_a_placement.mako"/>
 %endif

 %if attributes['gen_crt_b'] and not attributes['gen_crt_a']:
   <%include file="crt_b_placement.mako"/>
 %endif

</volume>
