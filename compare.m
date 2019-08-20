function a= compare(d,h)
b=(tanh((h-d)*1e12)+1)/2;
a=2-((tanh((b)*1e12)+1)/2)*2;