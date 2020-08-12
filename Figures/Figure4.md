# Figure 4 panels - October 2019
```R
## FIGURE4 PANEL B ABUNDANCE OF REPEATS
abundance_df = read.table(text = "
Assembly	LINE	LINEp	LTR	LTRp	Satellite	Satellitep	Simple_repeats	Simple_repeatsp	Low_complexity	Low_complexityp	Unknown	Unknownp	Total	Totalp
lycPyrIL	43000139	98.14427577	33335994	71.18048747	1012960	10.73582066	9066693	92.46812587	2296959	97.49831911	1413775	25.17717203	94042740	77.16605151
lycPyrSN1	42407197	96.79093449	35631863	76.08272841	947153	10.0383675	8175276	83.37686632	2050876	87.05290896	1413775	25.17717203	94646366	77.66135221
lycPyrSN2	43545882	99.38988922	40962836	87.46565753	1228750	13.02286333	8640208	88.11855005	2125889	90.23696292	1746137	31.09603129	101854789	83.5761686
lycPyrPB	43513457	99.31588183	46290202	98.84088483	9286689	98.42464423	9609444	98.00345918	2336184	99.16329074	5586383	99.48494338	120539720	98.90794592
lycPyrILPB	43228078	98.66452777	34462957	73.58682867	1300336	13.78156501	9538040	97.27523401	2450723	104.0250928	1558165	27.74853726	96456116	79.14632874
lycPyr2	43509695	99.30729538	46282990	98.82548545	9255441	98.09346341	9603513	97.94297092	2332154	98.99223056	5585744	99.47356377	120530085	98.90004
lycPyr4	43673424	99.68099333	46411817	99.10056255	9283945	98.39556204	9786408	99.80825498	2348562	99.68869594	5592762	99.59854362	121107744	99.37403367
lycPyr5	43795174	99.95887768	46797944	99.92503798	9426742	99.90899098	9844229	100.3979517	2387795	101.3540071	5642803	100.4896974	121903034	100.0266028
Final	43813191	100	46833051	100	9435329	100	9805209	100	2355896	100	5615305	100	121870613	100", sep = "\t", stringsAsFactors = F, header = T)

library(tidyr)
new_df = gather(data = abundance_df, key = Assembly)
names(new_df)[1] = "Repeat"
new_df$Assembly = rep(abundance_df$Assembly, 14)
new_df$value = as.numeric(new_df$value)

library(ggplot2)
library(RColorBrewer)

plot_df = new_df[grepl(pattern = "p$", x = new_df$Repeat),]
plot_df = plot_df[plot_df$Repeat != "Totalp",]
plot_df = plot_df[plot_df$Assembly != "Final",]
plot_df$Repeat = factor(x = plot_df$Repeat, levels = c("LINEp", "LTRp", "Satellitep", "Simple_repeatsp", "Low_complexityp", "Unknownp"), labels = c("LINE", "LTR", "Satellite", "Simple repeats", "Low complexity", "Unknown"))
plot_df$Assembly = factor(x = plot_df$Assembly, levels = c("lycPyrIL", "lycPyrSN1", "lycPyrSN2", "lycPyrPB", "lycPyrILPB", "lycPyr2", "lycPyr4", "lycPyr5"), labels = c("lycPyrIL", "lycPyrSN1", "lycPyrSN2", "lycPyrPB", "lycPyrILPB", "lycPyr2", "lycPyr4", "lycPyr5"))

#plot_df$Repeat = as.character(plot_df$Repeat)

Figure4_panelb = ggplot(data = plot_df) + geom_bar(aes(x = Repeat, y = value, fill = factor(Assembly)), stat = "identity", position = position_dodge2(width = 0.5)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1")
ggsave(filename = "Figure4_panelb_oct.pdf", plot = Figure4_panelb, device = "pdf", scale = 1, width = 15, height = 8, units = "cm", dpi = 300, limitsize = FALSE)


## FIGURE 4 PANEL C WHAT BREAKS THE GAPS

totalSummary = read.table(file = "Figure4c_data", sep = "\t", stringsAsFactors = F, header = T)

library(ggplot2)
library(RColorBrewer)

DISCARD = c("SINE", "tRNA")
totalSummary = totalSummary[!(totalSummary$Var1 %in% DISCARD),]
totalSummary$Var1 = factor(totalSummary$Var1, levels = c("Complex", "LINE", "LTR", "DNA", "rRNA", "Satellite", "Simple repeat", "Low complexity", "Unknown"))


plot = ggplot(data = totalSummary) + geom_bar(aes(x = Var1, y = Count, fill = Assembly), stat = "identity", position = position_dodge2(width = 0.5)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1")

ggsave(filename = "Figure4c_aug.pdf", plot = plot, device = "pdf", scale = 1, width = 15, height = 8, units = "cm", dpi = 300, limitsize = FALSE)


## FIGURE 4 PANEL D WHAT BREAKS THE GAPS PERCENTAGE OF GAPS
totalSummary = read.table(file = "Figure4c_data", sep = "\t", stringsAsFactors = F, header = T)
DISCARD = c("SINE", "tRNA")
totalSummary = totalSummary[!(totalSummary$Var1 %in% DISCARD),]
totalSummary$Var1 = factor(totalSummary$Var1, levels = c("Complex", "LINE", "LTR", "DNA", "rRNA", "Satellite", "Simple repeat", "Low complexity", "Unknown"))

gapIL = 14573
# not counting gaps smaller than 10
#gapIL = 7869
gapPB = 3422 * 2
gapSN1 = 21550
gapSN2 = 20131
gapDV = 6736

totalSummary$NGaps[totalSummary$Assembly == "lycPyr"] = gapIL
totalSummary$NGaps[totalSummary$Assembly == "lycPyr_SN1"] = gapSN1
totalSummary$NGaps[totalSummary$Assembly == "lycPyr_SN2"] = gapSN2
totalSummary$NGaps[totalSummary$Assembly == "lycPyr2"] = gapPB
totalSummary$NGaps[totalSummary$Assembly == "lycPyr2.1"] = gapDV

totalSummary$Freq = (totalSummary$Count/totalSummary$NGaps) * 100

plot = ggplot(data = totalSummary) + geom_bar(aes(x = Var1, y = Freq, fill = Assembly), stat = "identity", position = position_dodge2(width = 0.5)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1")

ggsave(filename = "Figure4d_aug.pdf", plot = plot, device = "pdf", scale = 1, width = 15, height = 8, units = "cm", dpi = 300, limitsize = FALSE)



## FIGURE4 GAP CONTENT
library(ggplot2)
library(RColorBrewer)

totalSummary = read.table("Figure4_gap_content_data", sep = "\t", stringsAsFactors = F, header = T)

DISCARD = c("DNA", "rRNA")

ggplot(data = totalSummary) + geom_bar(aes(x = Var1, y = Freq, fill = Assembly), stat = "identity", position = position_dodge2(width = 0.5)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1")


```
