{
  //ild TStyle
  TStyle* ildStyle = new TStyle("ildStyle", "ILD Style");

  //remove frame
  ildStyle -> SetFrameBorderMode(0);
  ildStyle -> SetCanvasBorderMode(0);
  ildStyle -> SetPadBorderMode(0);
  ildStyle -> SetTitleBorderSize(0);

  //use the primary color palette
  ildStyle -> SetPalette(1,0);

  //margins
  ildStyle -> SetPadTopMargin(0.13);
  ildStyle -> SetPadRightMargin(0.12);
  ildStyle -> SetPadLeftMargin(0.15);
  ildStyle -> SetPadBottomMargin(0.12);
  ildStyle -> SetCanvasColor(0);

  //fitting option
  ildStyle -> SetFuncColor(kRed);
  ildStyle -> SetFuncWidth(1);

  //axis label
  ildStyle -> SetLabelOffset(0.01,"xyz");
  ildStyle -> SetLabelColor(kBlack,"xyz");
  ildStyle -> SetLabelFont(62,"xyz");
  ildStyle -> SetLabelSize(0.05,"xyz");

  //title
  ildStyle -> SetTitleFont(62,"xyz");
  ildStyle -> SetTitleX(0.1);
  ildStyle -> SetTitleFillColor(0);
  ildStyle -> SetTitleFontSize(0.1);
  ildStyle -> SetTitleXOffset(1.1);
  ildStyle -> SetTitleYOffset(1.5);
  ildStyle -> SetTitleSize(0.05,"xyz");

  //stat
  ildStyle -> SetStatBorderSize(1);
  ildStyle -> SetStatH(0.2);
  ildStyle -> SetStatW(0.2);
  ildStyle -> SetStatX(0.9);
  ildStyle -> SetStatY(1);
  ildStyle -> SetStatFont(62);
  ildStyle -> SetStatColor(0);

  //legend
  ildStyle -> SetLegendBorderSize(1);
  ildStyle -> SetLegendFillColor();
  ildStyle -> SetLegendFont(62);
  ildStyle -> SetLegendTextSize(0.03);
  //histogram
  //ildStyle -> SetHistLineColor(kBlue);

  //the tick mark style
  ildStyle -> SetPadTickX(1);
  ildStyle -> SetPadTickY(1);
  ildStyle -> SetMarkerSize(0.7);

  ildStyle -> cd();
  gStyle -> ls();
}
