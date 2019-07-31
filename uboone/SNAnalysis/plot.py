# CAVEAT: adaptation of the SN Stream plot script for the trigger primitive analysis
# The SN stream role is played by the trigger primitives
# The EXT unbiased role is played by the gaushits (from the SN stream)

from ROOT import TFile, TH1F, TCanvas, gDirectory, gROOT, gPad, TLegend, TPad, TLine, TH1D, TLatex
import math

gROOT.ForceStyle()

### User-defined variables ###

# ROOT files with histograms
snfile = TFile("/home/jcrespo/MicroBooNE/TriggerDevelopment/analysis/SNRun19021_960filesPlaneY_TrigPrimRecoAna_slimmed_hist.root")
snname = "Trig. Prim."

extunbfile = snfile
extunbname = "gaushits"
extunbcolor = 2 # Red
extunbcolor_errors = 46 # Reddish

# Name of histograms to retrieve
hprimlistname = [ "hePrimSpectrum",
          "hgPrimSpectrum",
          "htotPrimSpectrum",
          "hePrimHitMult",
          "hgPrimHitMult",
          "htotPrimHitMult",
          "hePrimHitSpectrum",
          "hgPrimHitSpectrum",
          "htotPrimHitSpectrum",
          "hEventPrimHitMult",
          "hEventPrimHitSpectrum" ]

hlistname = [ "heSpectrum",
          "hgSpectrum",
          "htotSpectrum",
          "heHitMult",
          "hgHitMult",
          "htotHitMult",
          "heHitSpectrum",
          "hgHitSpectrum",
          "htotHitSpectrum",
          "hEventHitMult",
          "hEventHitSpectrum" ]

# Print number of entries in legend
printentries = False

# Fraction of canvas for the bottom pad (make it 0 if you do not want one)
frontierpad = 0.3
#frontierpad = 0
logx = False
#logx = True
logy = False
#logy = True

# Plot extrension
printext = ".png" # e.g. ".png"

# Print "MicroBooNE preliminary"
preliminary = False

# Event sizes in number of samples
snevtsize = 3200. # Only unique samples
extunbevtsize = 3200.

### End of user-defined variables ###

# Lists to hold histograms retrieved from files
snhlist= []
extunbhlist = []

# Normalization factor for the 3 primitive modes
# To do: change to single number 
normfactor = [0., 0., 0.]

# Retrieve objects from SN stream file
for rootdir in snfile.GetListOfKeys():
    if not rootdir.GetName().startswith("snmichel"):
        continue
    snfile.cd( rootdir.GetName() )
    for hname in hprimlistname:
        htemp = TH1F()
        gDirectory.GetObject( hname, htemp )
        htemp.SetNameTitle( "sn_" + htemp.GetName() + "_" + rootdir.GetName()[-1], htemp.GetTitle() + " tpm " + rootdir.GetName()[-1] )
        #htemp.SetNameTitle( "sn_" + htemp.GetName() + "_" + "2", htemp.GetTitle() + " plane " + "2" )
        htemp.Sumw2()
        snhlist.append( htemp )
        # Get normalization factor numerator
        if htemp.GetName().startswith("sn_hEventPrimHitMult"):
            # print "%s %i " % (htemp.GetName(), htemp.GetEntries())
            normfactor[int( rootdir.GetName()[-1] )] = htemp.GetEntries()*snevtsize
            #normfactor[2] = htemp.GetEntries()*snevtsize
            snexposure = htemp.GetEntries()*snevtsize*0.5e-6 # in seconds
            print "Trig. Prim. exposure: %i events = %.2f min " % (htemp.GetEntries(), snexposure/60.)

# Retrieve objects from EXT unbiased file
for rootdir in extunbfile.GetListOfKeys():
    if not rootdir.GetName().startswith("snmichel"):
        continue
    extunbfile.cd( rootdir.GetName() )
    for hname in hlistname:
        htemp = TH1F()
        gDirectory.GetObject( hname, htemp )
        # Get exposure (number of events) from hEventHitMult
        htemp.SetNameTitle( "extunb_" + htemp.GetName() + "_" + rootdir.GetName()[-1], htemp.GetTitle() + " tpm " + rootdir.GetName()[-1] )
        #htemp.SetNameTitle( "extunb_" + htemp.GetName() + "_" + "2", htemp.GetTitle() + " plane " + "2" )
        # Get normalization factor denominator
        if htemp.GetName().startswith("extunb_hEventHitMult"):
            # print "%s %i " % (htemp.GetName(), 2.*htemp.GetEntries())
            normfactor[int( rootdir.GetName()[-1] )] *= 1./(htemp.GetEntries()*extunbevtsize)
            #normfactor[2] *= 1./(htemp.GetEntries()*extunbevtsize)
            extunbexposure = htemp.GetEntries()*extunbevtsize*0.5e-6 # in seconds
            print "gaushit exposure: %i events = %.2f min " % (htemp.GetEntries(), extunbexposure/60.)
            # For hit multiplicity/event we have to correct for the different event size
            hnew = htemp.Clone()
            hnew.Reset()
            for ibin in xrange( htemp.GetNbinsX() ):
                # Trigger stream events are twice as long as SN stream ones, so they can have twice as many hits
                ibinnew = int( math.floor(ibin*snevtsize/extunbevtsize) )
                hnew.SetBinContent( ibinnew, hnew.GetBinContent(ibinnew) + htemp.GetBinContent(ibin) )
            htemp = hnew    
        htemp.Sumw2()
        extunbhlist.append( htemp )

# Turn off display of the EXT unbiased histograms
#normfactor[0] = 0.
#normfactor[1] = 0.
#normfactor[2] = 0.

# Lists to collect objects going out of scope
clist = []
extunbh_errorslist = []
leglist = []
hratiolist = []
line1list = []

# Loop for making plots
for snh in snhlist:

    c = TCanvas( "c_" + snh.GetName(), snh.GetTitle(), 800 if frontierpad else 1600, 800 )

    if frontierpad:
        ### Histograms plot ###
        pad1 = TPad( "pad1_" + snh.GetName(), snh.GetTitle(), 0., frontierpad, 1., 1. )
        pad1.SetBottomMargin(0) # Join with lower plot
        pad1.Draw()
        pad1.cd()

    extunbh = extunbhlist[snhlist.index(snh)]
    extunbh_integral = extunbh.Integral()
    # Apply normalization factor based on last character of name
    if extunbh.GetName().startswith("extunb_hEventHitMult"):
        # For the event hit multiplicity the event size is corrected differently
        extunbh.Scale( extunbevtsize/snevtsize*normfactor[int( snh.GetName()[-1] )] )
    else:
        extunbh.Scale( normfactor[int( snh.GetName()[-1] )] )
    extunbh.SetMarkerColor(extunbcolor)
    extunbh.SetLineColor(extunbcolor)
    extunbh_errors = extunbh.Clone( extunbh.GetName() + "_errors" )
    extunbh_errors.SetFillColor(extunbcolor_errors)
    extunbh_errors.Draw("E2") # Rectangle error boxes in the background
    #extunbh_errors.DrawNormalized("E2") # Rectangle error boxes in the background
    extunbh_errors.GetYaxis().SetRangeUser( 0.001, 
        1.1*extunbh_errors.GetMaximum() if ( extunbh_errors.GetMaximum() > snh.GetMaximum() ) else 1.1*snh.GetMaximum()  )
    #extunbh_errors.GetYaxis().SetLabelSize( (1. - frontierpad)*extunbh_errors.GetYaxis().GetLabelSize() )
    extunbh.Draw("hist same") # Central value line
    #extunbh.DrawNormalized("hist same") # Central value line

    snh.SetMarkerStyle(21)
    snh.SetMarkerSize(0.5)
    snh.Draw("E P0 same") # SN data
    #snh.DrawNormalized("E P0 X0 same") # SN data
    leg = TLegend(0.6, 0.7, 1.0, 1.0)
    leg.AddEntry( extunbh_errors, ("%s (%i entries)" % (extunbname, int( extunbh_integral ))) if printentries else "%s" % extunbname, "FLP" )
    leg.AddEntry( snh, ("%s (%i entries)" % (snname, int( snh.Integral() ))) if printentries else "%s" % snname, "LEP" )
    if logx: leg.SetFillStyle(0) # Make it transparent to mitigate the change of the distribution
    leg.Draw()
    if logx: gPad.SetLogx()
    if logy: 
        extunbh_errors.GetYaxis().SetRangeUser( 0.9, 1.1*extunbh_errors.GetMaximum() )
        gPad.SetLogy()

    if preliminary:
        tx = TLatex()
        tx.SetTextSize(0.04)
        tx.SetTextAlign(11) # Bottom left adjusted
        tx.DrawTextNDC( 0.15, 0.95, "MicroBooNE Preliminary");

    gPad.Modified()
    gPad.Update()
    print "%s: rate %.2f +/- %.2f" % ( snh.GetName(), snh.Integral()/snexposure, math.sqrt(snh.Integral())/snexposure )
    print "\t trigger: rate %.2f +/- %.2f" % ( extunbh_integral/extunbexposure, math.sqrt(extunbh_integral)/extunbexposure )

    ### Ratio plot ###
    if frontierpad:
        c.cd()
        pad2 = TPad( "pad2_" + snh.GetName(), snh.GetTitle(), 0., 0.05, 1., frontierpad )
        pad2.SetTopMargin(0) # Join with upper plot
        pad2.SetBottomMargin(0.3)
        pad2.Draw()
        pad2.cd()

        hratio = snh.Clone()
        hratio.Sumw2()
        hratio.Divide(extunbh)
            
        hratio.GetYaxis().SetTitle("%s/%s" % (snname, extunbname))
        hratio.GetYaxis().SetRangeUser(0., 1.999)
        hratio.GetYaxis().SetNdivisions(505)
        hratio.GetYaxis().SetLabelSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetYaxis().GetLabelSize() )
        hratio.GetYaxis().SetTitleSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetYaxis().GetTitleSize() )
        hratio.GetYaxis().SetTitleOffset( 1.10*(frontierpad/(1. - frontierpad)*extunbh_errors.GetYaxis().GetTitleOffset()) )
        hratio.GetYaxis().CenterTitle()
        hratio.GetXaxis().SetLabelSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetXaxis().GetLabelSize() )
        hratio.GetXaxis().SetTitleSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetXaxis().GetTitleSize() )
        hratio.SetMarkerStyle(21)
        hratio.SetMarkerSize(0.5)
        hratio.Draw("E P0")
        line1 = TLine( extunbh_errors.GetXaxis().GetXmin(), 1., extunbh_errors.GetXaxis().GetXmax(), 1. )
        line1.SetLineWidth(1)
        line1.SetLineColor(extunbcolor)
        line1.Draw()
        hratio.Draw("E P0 same") # Draw again over line
        if logx: gPad.SetLogx()        
        # Keep the objects in memory
        line1list.append(line1)
        hratiolist.append(hratio)

    leglist.append(leg)
    extunbh_errorslist.append(extunbh_errors)
    c.Print( c.GetName() + ("_logx" if logx else "") + ("_logy" if logy else "") + ("_ratio" if frontierpad else "") + printext ) 
    #c.Print( c.GetName() + ("_logx" if logx else "") + ("_bigratio" if frontierpad else "") + printext)
    #c.Print( "norm_" + c.GetName() + printext )
    clist.append(c)

raw_input("Press the <ENTER> key to continue...")
