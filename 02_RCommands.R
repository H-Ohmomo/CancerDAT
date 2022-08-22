###########################################################################
## 本Rcommandsファイルには、本章で使用するRコマンドのみ含まれております．##
## Webサイトからダウンロード等操作の詳細は本書をご確認ください．			 	 ##
###########################################################################

# MacのアプリケーションからRStudioを起動し、各種パッケージをインストールする．
install.packages("BiocManager")
library(BiocManager)
BiocManager::install(c("minfi", "minfiData", "ggplot2"), force = TRUE)

# インストールした各種パッケージを呼び出す．
library(minfi)
library(minfiData)
library(ggplot2)

# 解析データを、「第８章がんのエピゲノム解析、8.2  解析環境構築編、④  析に使用する公 データをダウンロードする」を参考に準備する．
# サンプルシートを、「第８章がんのエピゲノム解析、8.2  解析環境構築編、⑤  析データのサンプルシートを作成する」を参考に準備する．

# 解析の前準備を行う．
setwd("~/Desktop/DNAm_analysis")
data <- read.csv("E-GEOD-61441_SS2.csv", header=T)
data
udata <- unique(data)
udata
write.csv(udata, row.names=FALSE, "E-GEOD-61441.raw.all/E-GEOD-61441_SS3.csv")

# 作成したサンプルシートをテキストエディタで開いて、最後の行に改行を入れる．
# HM450のマニフェストファイルを、「第８章がんのエピゲノム解析、8.2  解析環境構築編、 析に使用するマニフェストファイルをダウンロードする」を参考に準備する．

# 8.3 DNAメチル化率算出編
# ① Rを起動し  析に使用するパッケージを呼び出す
setwd("~/Desktop/DNAm_analysis")
install.packages("BiocManager")　 #すでに実施済みの場合不要
library(BiocManager)　#すでに実施済みの場合不要
BiocManager::install(c("minfi", "minfiData", "ggplot2"), force = TRUE) 　#すでに実施済みの場合不要
library(minfi)　#すでに実施済みの場合不要
library(minfiData)　#すでに実施済みの場合不要
library(ggplot2)　#すでに実施済みの場合不要

# 解析を行うディレクトリを指定し、サンプルシートを読み込む
basedir <- "E-GEOD-61441.raw.all"
targets <- read.metharray.sheet(basedir)
targets

# HM450データ(idatファイル)を読みだし、RGsetにおさめる．
RGset <- read.metharray.exp(targets=targets, force=TRUE)
RGset
save(RGset,file="02_E-GEOD-61441_RGset.RData")

# ② データの正 化およびβ値 メチル化率 を算出する
mset.noob <- preprocessNoob(RGset)
save(mset.noob, file = "03_E-GEOD-61441_mset.RData")
beta.noob <- getBeta(mset.noob)
save(beta.noob, file = "04_E-GEOD-61441_beta.RData")

# ③ コントロールプローブの蛍光強度を取得，Principal Conponent Analysis，はずれ値サンプルがないかを確認
## コントロールプローブの情報を取り出し、まとめる．
ctrls <- getProbeInfo(RGset, type = "Control")
ctrls <- ctrls[ctrls$Address %in% featureNames(RGset),]
ctrls

ctrl_r <- getRed(RGset)[ctrls$Address[!(ctrls$Type %in% "NEGATIVE")],]
ctrl_g <- getGreen(RGset)[ctrls$Address[!(ctrls$Type %in% "NEGATIVE")],]
ctrl_nneg = rbind(ctrl_r, ctrl_g)
save(ctrl_nneg, file = "05_E-GEOD-61441_ctrlnneg.RData")

pctrl <- ctrls$ExtendedType[!(ctrls$Type %in% "NEGATIVE")]
r_pctrl <- paste(pctrl,"R",sep="")
g_pctrl <- paste(pctrl,"G",sep="")
ct <- c(r_pctrl, g_pctrl)

neg <- t(ctrl_nneg)
colnames(neg) <- ct

SampleID <- targets[,1]
neg2 <- cbind(SampleID,neg)
write.table(neg2,"05_E-GEOD-61441_ctrlprobe.txt",row.name=F,quote=F)

## コントロールプローブの情報を使ってPCAを行う．
HM450_ctrlpca <- prcomp(neg)
save(HM450_ctrlpca, file = "06_E-GEOD-61441_ctrlpca.RData")

options(scipen=10)
p <- HM450_ctrlpca$x
cpca <- cbind(SampleID,p)
write.table(cpca,"07_E-GEOD-61441_ctrlprobepca_sorted.txt",row.name=F,quote=F)

## コントロールプローブの情報を使って作図する．
## 作図1. Sample_Groupで色分けしたPCAを描画する．
fn.sample.pca <- paste("07_E-GEOD-61441_ctrlp_pca.png",sep="")

ggplot.pca.point <- function(fname, pca, pc1, pc2, Sample_Group) {
	pca.df <- data.frame(
		pc1 = pca$x[,pc1],
		pc2 = pca$x[,pc2],
		sample_group = Sample_Group,
		label = attr(pca$x, "dimnames")[[1]])
		
		xlabel <- paste("PC", pc1, " (", summary(pca)$importance[2,][pc1]*100, "%)", sep="")
		ylabel <- paste("PC", pc2, " (", summary(pca)$importance[2,][pc2]*100, "%)", sep="")

		g <- ggplot(pca.df, aes(x=pc1, y=pc2, colour = Sample_Group))
		g <- g + geom_point()
		g <- g + labs(x = xlabel, y = ylabel)
		g <- g + xlim(-90000,90000)+ylim(-50000,50000)
		plot(g)
		ggsave(filename=fname, width=8, height=8, dpi=300)
}

ggplot.pca.point(fn.sample.pca, HM450_ctrlpca, 1, 2, factor(targets$Sample_Group))


 # ④ 性別推定，フィルタリングの設定を行い，常染色体のみのEWAS用のデータセットを整形する．
Gmset.noob <- mapToGenome(mset.noob)
sex <- getSex(Gmset.noob)
sex
save(sex,file = "08_E-GEOD-61441_predsex.RData")
est_sex_res <- data.frame(targets[, c("SampleID", "Sample_Sex")], sex, row.names = NULL)
est_sex_res$Sample_Sex <- as.factor(est_sex_res$Sample_Sex)
est_sex_res$predictedSex <- factor(est_sex_res$predictedSex, levels = c("F", "M"), labels = levels(est_sex_res$Sample_Sex))
est_sex_res
identical(targets$SampleID, rownames(sex))

# targetsにおさめられているSampleIDと sexにおさめられているサンプル名が一致しているか作図して確認する．
pdf("predictedSex.pdf")
plot(
	x = est_sex_res$xMed, xlab = "xMed",
	y = est_sex_res$yMed, ylab = "yMed",
	pch = 19,
	col = c("#FF6A6A", "#009ACD")[est_sex_res$Sample_Sex]
	)

abline(a = -2, b = 1, lty = 2, col = "#919191")
abline(a = -3, b = 1, lty = 1, col = "#919191")

text(
	x = est_sex_res$xMed, xlab = "xMed",
	y = est_sex_res$yMed, ylab = "yMed",
	labels = sex$predictedSex,
	col = c("#FF6A6A", "#009ACD")[est_sex_res$predictedSex],
	pos = 4
	)

legend(
	"bottomleft",
	title = "Sample_Sex",
	legend = c(levels(est_sex_res$Sample_Sex), paste("cutoff =", c(-2, -3))),
	bg = "#FFFFFF",
	pch = c(19, 19, NA, NA),
	lty = c(0, 0, 2, 1),
	col = c("#FF6A6A", "#009ACD", "#919191", "#919191"),
	inset = .01
	)
dev.off()

# detection p valueによるフィルタリング閾値を設定して、外れ値を本解析に含めないように処理する．
pval <- detectionP(RGset,m+U)
save(pval,file="09_E-GEOD-61441_detectionP.RData")
p2 <- pval[rownames(beta.noob),]
colnames(p2) <- colnames(beta.noob)
thres = 1E-16
b = ifelse(p2<thres,beta.noob,NA)

## サンプルごとのコールレートを算出する．
sample.call=colSums(!is.na(b))/nrow(b)

## CpGごとのコールレートを算出する．
marker.call=rowSums(!is.na(b))/ncol(b)

## コールレートの閾値は一般的な0.95(95%)とする．
callrate.thres = 0.95

## サンプルのコールレートでフィルタリングする．
samples=names(sample.call[sample.call>callrate.thres])

## CpGごとのコールレートでフィルタリングする．
markers=names(marker.call[marker.call>callrate.thres])
save(sample.call,marker.call,file="10_E-GEOD-61441_callrate.RData")

## フィルタリングを行ったβ値をb2に収めて保存する．
b2 <- beta.noob[markers,samples]
save(b2,file="11_E-GEOD-61441_betaf.RData")

## b2に収められた結果を確認する．
b2


# マニフェストファイルの情報をもとに常染色体のみのβ値のデータセットに整形する．
path.manifest <- "HumanMethylation450_15017482_v1-2.csv"
anno <- read.csv(path.manifest, as.is=TRUE, skip = 7)
anno2 = anno[,c('Infinium_Design_Type', 'Color_Channel', 'CHR', 'MAPINFO', 'Strand', 'Name')]
cas = anno[substr(anno2$Name, 1,3)=='ch.' & !(anno2$CHR %in% c('X','Y')),]
cgs = anno[substr(anno2$Name, 1,2)=='cg' & !(anno2$CHR %in% c('X','Y')),]
auto = c(cgs$Name, cas$Name)
auto = as.matrix(auto)
auto2 <- c(auto)
beta <- b2[row.names(b2) %in% auto2,]
sample <- subset(targets, select=c(SampleID))
colnames(beta) <- sample$SampleID
probeID <- rownames(beta)
p <- cbind(probeID,beta)
write.table(p,"12_E-GEOD-61441_beta.csv",row.name=F,quote=F)

## メチル化率がNAとなっているCpGサイトを除去し、常染色体のβ値を使ってPCAを作図する．
betat <- t(beta)
betat.naomit <- na.omit(betat)
pca.betat.naomit.auto <- prcomp(betat.naomit)
head(pca.betat.naomit.auto)
save(pca.betat.naomit.auto,file="13_E-GEOD-61441_betafpcaauto.RData")
autop <- pca.betat.naomit.auto$x
cautopca <- cbind(SampleID, autop)
write.table(cautopca,"14_E-GEOD-61441_cautopca.txt",row.name=F,quote=F)

## 作図1. Sample_Groupで色分けしたPCAを描画する．
fn.sample.pca <- paste("E-GEOD-61441_Group_autob_pca.png",sep="")

ggplot.pca.point <- function(fname, pca, pc1, pc2, Sample_Group) {
	pca.df <- data.frame(
		pc1 = pca$x[,pc1],
		pc2 = pca$x[,pc2],
		sample_group = Sample_Group,
		label = attr(pca$x, "dimnames")[[1]])
		xlabel <- paste("PC", pc1, " (", summary(pca)$importance[2,][pc1]*100, "%)", sep="")
		ylabel <- paste("PC", pc2, " (", summary(pca)$importance[2,][pc2]*100, "%)", sep="")
		g <- ggplot(pca.df, aes(x=pc1, y=pc2, colour=Sample_Group))
		g <- g + geom_point()
		g <- g + labs(x = xlabel, y = ylabel)
		g <- g + xlim(-100,100)+ylim(-50,50)
		plot(g)
		ggsave(filename=fname, width=8, height=8, dpi=300)
}

ggplot.pca.point(fn.sample.pca, pca.betat.naomit.auto, 1, 2, factor(targets$Sample_Group))


## 作図2. Sample_Sexで色分けしたPCAを描画する．
fn.sample.pca2 <- paste("E-GEOD-61441_Sex_autob_pca.png",sep="")

ggplot.pca.point <- function(fname, pca, pc1, pc2, Sample_Sex) {
	pca.df <- data.frame(
		pc1 = pca$x[,pc1],
		pc2 = pca$x[,pc2],
		sample_sex = Sample_Sex,
		label = attr(pca$x, "dimnames")[[1]])
		xlabel <- paste("PC", pc1, " (", summary(pca)$importance[2,][pc1]*100, "%)", sep="")
		ylabel <- paste("PC", pc2, " (", summary(pca)$importance[2,][pc2]*100, "%)", sep="")
		g <- ggplot(pca.df, aes(x=pc1, y=pc2, colour=Sample_Sex))
		g <- g + geom_point()
		g <- g + labs(x = xlabel, y = ylabel)
		g <- g + xlim(-100,100)+ylim(-50,50)
		plot(g)
		ggsave(filename=fname, width=8, height=8, dpi=300)
}

ggplot.pca.point(fn.sample.pca2, pca.betat.naomit.auto, 1, 2, factor(targets$Sample_Sex))


## 作図3. Sample_PairedIDで色分けしたPCAを描画する．
fn.sample.pca3 <- paste("E-GEOD-61441_PairedID_autob_pca.png",sep="")

ggplot.pca.point <- function(fname, pca, pc1, pc2, Sample_PairedID) {
	pca.df <- data.frame(
		pc1 = pca$x[,pc1],
		pc2 = pca$x[,pc2],
		sample_pairedid = Sample_PairedID,
		label = attr(pca$x, "dimnames")[[1]])
		xlabel <- paste("PC", pc1, " (", summary(pca)$importance[2,][pc1]*100, "%)", sep="")
		ylabel <- paste("PC", pc2, " (", summary(pca)$importance[2,][pc2]*100, "%)", sep="")
		g <- ggplot(pca.df, aes(x=pc1, y=pc2, colour=Sample_PairedID))
		g <- g + geom_point()
		g <- g + labs(x = xlabel, y = ylabel)
		g <- g + xlim(-100,100)+ylim(-50,50)
		plot(g)
		ggsave(filename=fname, width=8, height=8, dpi=300)
}

ggplot.pca.point(fn.sample.pca3, pca.betat.naomit.auto, 1, 2, factor(targets$Sample_PairedID))


# 8.4 エピゲノム関連解析（EWAS）編
# ① EWASに使用するパッケージをインストールまたは呼び出す
setwd("~/Desktop/DNAm_analysis")　#8.3 DNAメチル化率算出編から通しで解析している場合は不要です．
BiocManager::install(c("data.table", "qqman", "dplyr", "ggplot2"), force = TRUE)
library(data.table)
library(qqman)
library(dplyr)
library(ggplot2)

#② 先の解析で作成したβ値データのファイルをEWAS析用に整える
beta=as.data.frame(fread("12_E-GEOD-61441_beta.csv",header=T),stringsAsFactors=F)
cov=read.csv("E-GEOD-61441.raw.all/E-GEOD-61441_SS3.csv",header=T,stringsAsFactors=F)
id_beta = colnames(beta)[2:ncol(beta)]
id_cov = cov$SampleID
ids = intersect(id_beta, id_cov)
beta2 = beta[,c("probeID",ids)]
rownames(cov) = cov$SampleID
cov2 = cov[ids,]
ftable(cov2[,c("Sample_Group","Sample_Sex")])

##横軸にPC1を、縦軸にPC2でPCAを描画する．
tmp1 = na.omit(beta2[,c(2:ncol(beta2))])
sd = apply(tmp1, 1, sd)
tmp2 = tmp1[sd 0,]
pca = prcomp(t(tmp2),scale=T)
tmp3 = as.data.frame(t(summary(pca)$importance))[,c("Proportion of Variance","Cumulative Proportion")]
tmp4 = as.data.frame(pca$x)
tmp4$SampleID = rownames(tmp4)
tmp5 = merge(tmp4,cov2,by="SampleID",all=F)
rownames(tmp5)=tmp5$SampleID
cov3=tmp5[ids,]

## PC値01～92をサンプルシートに挿入して 出力する．
write.table(cov3, "pcaRes.tsv", row.names=F, quote=F, sep="\t")

## Proportion of variance（寄与率）とCumulative proportion（累積寄与率）を出力する．
tmp3$PC = rownames(tmp3)
write.table(tmp3, "pcaProportion.tsv", row.names=F, quote=F, sep="\t")

## PC1とPC2の期待値をそれぞれexpPC1とexpPC2に代入する．
expPC1 = round(summary(pca)$importance[2,"PC1"],digits=4)
expPC2 = round(summary(pca)$importance[2,"PC2"],digits=4)

g = ggplot(tmp5, aes(x = PC1, y = PC2, colour = Sample_Group)) + geom_point() + theme_bw() + xlab(paste0("PC1(", 100*expPC1, "%)")) + ylab(paste0("PC2(", 100*expPC2, "%)"))

ggsave("pca1-2.pdf", g, height=5, width=6)

## 凡例を示さない図を描画するには、 theme(legend.position="NONE")を加える．
g2 = ggplot(tmp5, aes(x=PC1, y=PC2, colour = Sample_Group)) + geom_point() + theme_bw() + xlab(paste0("PC1 (", 100*expPC1, "%)")) + ylab(paste0("PC2(", 100*expPC2, "%)")) + theme(legend.position = "NONE")

ggsave("pca1-2_woLegend.pdf", g2, height=5, width=5)

#③ Sample_GroupとSample_Sexを補正項に入れてEWASを行う．
## EWASの式を作成する．
model = as.formula("m ~ Sample_Group+Sample_Sex")
ewas = function(x){
	probeID = as.character(x[1])
	m = as.numeric(x[2:length(x)])
	cr = (100*length(m[!is.na(m)]))/length(m)
	res1 = c(probeID, cr)
	tmp = cov3
	tmp$m = m
	tryCatch({tmp1 = lm(model, data=tmp, na.action=na.omit)}, error = function(e){tmp1 <<- NULL})
	if(is.null(tmp1)){
		res2 = c(rep(NA, 3))
	} else if(!is.na(tmp1$coefficients[2])){
	tmp2 = summary(tmp1)$coefficients
	coe = tmp2[2, 1]
	pval = tmp2[2, 4]
	rsq = summary(tmp1)$r.squared
	res2 = c(coe, pval, rsq)
	} else {
	res2 = c(rep(NA, 3))
	}
	return(c(res1, res2))
}

## EWASを行って結果を出力する．
ewasRes = as.data.frame(t(apply(beta2, 1, ewas)), stringsAsFactors=F)
colnames(ewasRes) = c("probeID", "cr", "coef", "pval", "rsq")
write.table(ewasRes, "ewasRes.tsv", row.names=F, quote=F, sep="\t")

## EWAS結果に、マニフェストファイルから遺伝子名などをつける．
manif = read.csv("HumanMethylation450_15017482_v1-2.csv", as.is=TRUE, skip = 7)
annot = data.frame(probeID=manif$IlmnID, chr=manif$CHR, pos=manif$MAPINFO, genename=manif$UCSC_RefGene_Name, stringsAsFactors=F)
ewasRes2 <- inner_join(ewasRes, annot, by = "probeID")
write.table(ewasRes2, "ewasRes2merged.tsv", row.names=F, quote=F, sep="\t")

## EWAS結果をもとにマンハッタンプロットを描画する．
ewasRes2$pval=as.numeric(ewasRes2$pval)
ewasRes2$pos=as.numeric(ewasRes2$pos)
ewasRes2$chr=as.numeric(ewasRes2$chr)
png("manhattan_ewasRes2.png")
manhattan(x=ewasRes2, snp = "probeID", chr = "chr", bp = "pos", p = "pval", suggestiveline = -log10(1e-05), genomewideline = -log10(1.06e-07))
dev.off()

## EWAS結果をもとにquantile‒quantile（QQ）プロットを描画する．
png("qqplot_ewasRes2.png")
qq(ewasRes2$pval)
dev.off()

## genomic inflation factor（λ（ラムダ）値）を算出する．
pvalue <- as.data.frame(as.numeric(ewasRes2$pval))
colnames(pvalue) <- "pval"

for(i in 1:ncol(pvalue)){
	pvalue[,i] <- sort(pvalue[,i],na.last=T)
}

pvalue$exp <- (1:nrow(pvalue) - 0.5) / nrow(pvalue)
pval.log <- -log10(pvalue)

Quant = function(List, thres){
	x = sort(List, decreasing = F, na.last = T)
	N = length(List)
	r = floor(N*thres-qnorm(0.975)*(N*thres*(1-thres))^0.5)
	s = ceiling(1+N*thres+qnorm(0.975)*(N*thres*(1-thres))^0.5)
	rep = c(quantile(x, probs = thres, na.rm = T), x[r], x[s])
	qchisq(rep, df=1, lower.tail = F)/qchisq(0.5, df=1)
}

lambda <- apply(pvalue, 2, function(x){Quant(x, 0.5)})
rownames(lambda) <- c("lambda", "confL", "confR")
lambda2 = paste0(round(lambda[1, 1],digits=3)," (",round(lambda[3, 1],digits=3), "-", round(lambda[2, 1],digits=3),")")

write.table(lambda2, row.names=FALSE, "lambda2.txt")

#④ PC1～20を補正項に追加してEWASを行う．
## EWASの式を作成する．
model2 = as.formula("m ~ Sample_Group + Sample_Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")
ewas2 = function(x){
	probeID = as.character(x[1])
	m = as.numeric(x[2:length(x)])
	cr = (100*length(m[!is.na(m)]))/length(m)
	res1 = c(probeID, cr)
	tmp = cov3
	tmp$m = m
	tryCatch({tmp1 = lm(model2, data=tmp, na.action=na.omit)}, error = function(e){tmp1 <<- NULL})
	if(is.null(tmp1)){
		res2 = c(rep(NA, 3))
	} else if(!is.na(tmp1$coefficients[2])){
	tmp2 = summary(tmp1)$coefficients
	coe = tmp2[2, 1]
	pval = tmp2[2, 4]
	rsq = summary(tmp1)$r.squared
	res2 = c(coe, pval, rsq)
	} else {
	res2 = c(rep(NA, 3))
	}
	return(c(res1, res2))
}

## EWASを行って結果を出力する．
ewasRes3 = as.data.frame(t(apply(beta2, 1, ewas2)), stringsAsFactors = F)
colnames(ewasRes3) = c("probeID", "cr", "coef", "pval", "rsq")
write.table(ewasRes3, "19_ewasRes3.tsv", row.names = F, quote = F, sep = "\t")

## EWAS結果に、マニフェストファイルから遺伝子名などをつける．
ewasRes4 <- inner_join(ewasRes3, annot, by = "probeID")
write.table(ewasRes4, "20_ewasRes4merged.tsv", row.names=F, quote=F, sep="\t")

## EWAS結果をもとにマンハッタンプロットを描画する．
ewasRes4$pval=as.numeric(ewasRes4$pval)
ewasRes4$pos=as.numeric(ewasRes4$pos)
ewasRes4$chr=as.numeric(ewasRes4$chr)
png("manhattan_ewasRes4.png")
manhattan(x=ewasRes4, snp = "probeID", chr = "chr", bp = "pos", p = "pval", suggestiveline = -log10(1e-05), genomewideline = -log10(1.06e-07))
dev.off()

## EWAS結果をもとにquantile‒quantile（QQ）プロットを描画する．
png("qqplot_ewasRes4.png")
qq(ewasRes4$pval)
dev.off()

## genomic inflation factor（λ（ラムダ）値）を算出する．
pvalue2 <- as.data.frame(as.numeric(ewasRes4$pval))
colnames(pvalue2) <- "pval"
for(i in 1:ncol(pvalue2)){
	pvalue2[,i] <- sort(pvalue2[,i],na.last=T)
}

pvalue2$exp <- (1:nrow(pvalue2) - 0.5) / nrow(pvalue2)
pvalue2.log <- -log10(pvalue2)

Quant2 = function(List, thres){
	x2 = sort(List, decreasing = F, na.last = T)
	N2 = length(List)
	r2 = floor(N2*thres-qnorm(0.975)*(N2*thres*(1-thres))^0.5)
	s2 = ceiling(1+N2*thres+qnorm(0.975)*(N2*thres*(1-thres))^0.5)
	rep2 = c(quantile(x2, probs = thres, na.rm = T), x2[r2], x2[s2])
	qchisq(rep2, df=1, lower.tail = F)/qchisq(0.5, df=1)
}

lambda3 <- apply(pvalue2, 2, function(x2){Quant2(x2, 0.5)})
rownames(lambda3) <- c("lambda", "confL", "confR")
lambda4 = paste0(round(lambda3[1, 1],digits=3)," (",round(lambda3[3, 1],digits=3), "-", round(lambda3[2, 1],digits=3),")")

write.table(lambda4, row.names=FALSE, "lambda4.txt")

