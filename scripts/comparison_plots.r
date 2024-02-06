library(ggplot2)

iupred_long_diso_table = read.table("consolidated_iupred-long.tsv", header=T)
iupred_short_diso_table = read.table("consolidated_iupred-short.tsv", header=T)
disopred_diso_table = read.table("consolidated_disopred.tsv", header=T)

iupred_long_plot <- ggplot(iupred_long_diso_table, aes(x=effector_type, y=disorder_fraction, color=collection_type)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA,linewidth = 0.55,size=.3)+
  geom_point(shape = 1, size=2/3,position = position_jitterdodge(jitter.width = .05,dodge.width = 0.9),alpha=0.8) +
  scale_color_manual(values=c("#e5be97", "#97bee5")) +
  xlab(label = "") +
  ylab(label = "Disorder fraction") +
  theme_classic() +
  theme(legend.position="top") +
  theme(legend.title=element_blank())

iupred_short_plot <- ggplot(iupred_short_diso_table, aes(x=effector_type, y=disorder_fraction, color=collection_type)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA,linewidth = 0.55,size=.3)+
  geom_point(shape = 1, size=2/3,position = position_jitterdodge(jitter.width = .05,dodge.width = 0.9),alpha=0.8) +
  scale_color_manual(values=c("#e5be97", "#97bee5")) +
  xlab(label = "") +
  ylab(label = "Disorder fraction") +
  theme_classic() +
  theme(legend.position="top") +
  theme(legend.title=element_blank())

disopred_plot <- ggplot(disopred_diso_table, aes(x=effector_type, y=disorder_fraction, color=collection_type)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA,linewidth = 0.55,size=.3)+
  geom_point(shape = 1, size=2/3,position = position_jitterdodge(jitter.width = .05,dodge.width = 0.9),alpha=0.8) +
  scale_color_manual(values=c("#e5be97", "#97bee5")) +
  xlab(label = "") +
  ylab(label = "Disorder fraction") +
  theme_classic() +
  theme(legend.position="top") +
  theme(legend.title=element_blank())

print(iupred_long_plot)
print(iupred_short_plot)
print(disopred_plot)

ggsave(filename= "comparison_iupred-long.pdf", plot=iupred_long_plot, dpi=300)
ggsave(filename= "comparison_iupred-long.svg", plot=iupred_long_plot, dpi=300)

ggsave(filename= "comparison_iupred-short.pdf", plot=iupred_short_plot, dpi=300)
ggsave(filename= "comparison_iupred-short.svg", plot=iupred_short_plot, dpi=300)

ggsave(filename= "comparison_disopred.pdf", plot=disopred_plot, dpi=300)
ggsave(filename= "comparison_disopred.svg", plot=disopred_plot, dpi=300)
