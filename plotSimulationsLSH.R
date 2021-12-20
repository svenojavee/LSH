library(data.table)
library(ggplot2)

mseDat <- fread("/Users/admin/Documents/plotLSH/mseDat.txt")
mseDat[variable=="origH2",variable:="Original expression"]
mseDat[variable=="newH2",variable:="New expression"]


p_0 <- ggplot(dat = mseDat[Imprecision == 0,], aes(x = K, y = mse, color = variable)) + geom_line() + facet_wrap(. ~ h2, scales = "free",ncol=3)+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 10)) +
    theme(axis.title = element_text(color = "gray20", size = 12)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 8, color = "gray10")) +
    xlab("Prevalence") + ylab("MSE") + scale_x_continuous(trans='log10')


p_01 <- ggplot(dat = mseDat[Imprecision == 0.1,], aes(x = K, y = mse, color = variable)) + geom_line() + facet_wrap(. ~ h2, scales = "free",ncol=3)+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 10)) +
    theme(axis.title = element_text(color = "gray20", size = 12)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 8, color = "gray10")) +
    xlab("Prevalence") + ylab("MSE") + scale_x_continuous(trans='log10')


p_02 <- ggplot(dat = mseDat[Imprecision == 0.2,], aes(x = K, y = mse, color = variable)) + geom_line() + facet_wrap(. ~ h2, scales = "free",ncol=3)+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, 'cm')) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 12), legend.key.size = unit(0.5, 'cm'), legend.spacing.x = unit(0.2, 'cm')) +
    theme(title = element_text(color = "gray20", size = 10)) +
    theme(axis.title = element_text(color = "gray20", size = 12)) +
    theme(axis.text = element_text(color = "gray40", size = 10)) +
    theme(legend.position = "top") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(size = 8, color = "gray10")) +
    xlab("Prevalence") + ylab("MSE") + scale_x_continuous(trans='log10') 

ggsave(p_0,filename="/Users/admin/Documents/plotLSH/p0.pdf",width=8,height=5)
ggsave(p_01,filename="/Users/admin/Documents/plotLSH/p_01.pdf",width=8,height=5)
ggsave(p_02,filename="/Users/admin/Documents/plotLSH/p_02.pdf",width=8,height=5)
