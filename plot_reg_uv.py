qplt.pcmeshclf(regu_low,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, u lower tropos')
plt.savefig('aus_ulow_reg.pdf')
qplt.pcmeshclf(coru_low,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, u lower tropos')
plt.savefig('aus_ulow_cor.pdf')
qplt.pcmeshclf(regu_mid,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, u mider tropos')
plt.savefig('aus_umid_reg.pdf')
qplt.pcmeshclf(coru_mid,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, u mider tropos')
plt.savefig('aus_umid_cor.pdf')
qplt.pcmeshclf(regu_up,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, u uper tropos')
plt.savefig('aus_uup_reg.pdf')
qplt.pcmeshclf(coru_up,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, u uper tropos')
plt.savefig('aus_uup_cor.pdf')

qplt.pcmeshclf(regv_low,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, v lower tropos')
plt.savefig('aus_vlow_reg.pdf')
qplt.pcmeshclf(corv_low,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, v lower tropos')
plt.savefig('aus_vlow_cor.pdf')
qplt.pcmeshclf(regv_mid,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, v mider tropos')
plt.savefig('aus_vmid_reg.pdf')
qplt.pcmeshclf(corv_mid,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, v mider tropos')
plt.savefig('aus_vmid_cor.pdf')
qplt.pcmeshclf(regv_up,cmap=mc.jetwhite())
plt.title('Regression between Tsfc Aus, v uper tropos')
plt.savefig('aus_vup_reg.pdf')
qplt.pcmeshclf(corv_up,vmin=-0.7,vmax=0.7,cmap=mc.jetwhite())
plt.title('Correlation between Tsfc Aus, v uper tropos')
plt.savefig('aus_vup_cor.pdf')

