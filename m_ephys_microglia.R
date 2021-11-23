# SETUP ----
library(dplyr)
library(ggplot2)
library(xlsx)
library(stats)
library(tidyverse)
library(fastDummies)
library(here)


## LOAD DATA ----
m_ephys <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/ephys/processed/mouse_ephys_features_gouwens.csv")

# get ephys session id from ephys data to match the metadata
m_ephys$ephys2 <- substr(m_ephys$ephys_session_id, nchar(m_ephys$ephys_session_id) - 21 + 1, nchar(m_ephys$ephys_session_id))
m_ephys$ephys3 <- substr(m_ephys$ephys2, 1, 9)
m_ephys$ephys_session_id <- as.integer(m_ephys$ephys3)

meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_patchseq_2020/mouse/20200711_patchseq_metadata_mouse.csv")
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)
meta$transcriptomics_sample_id = sub("-",".", meta$transcriptomics_sample_id)

micro <- read.xlsx(here("m_microglia_ranked.xlsx", sheetName = 'm_microglia_ranked'))

# merge
temp <- inner_join(meta, micro, by = "transcriptomics_sample_id")
mm <- dplyr::select(temp, cell_specimen_id, ephys_session_id, transcriptomics_sample_id, Micro.PVM)

final <- inner_join(mm, m_ephys, by = "ephys_session_id") # got 3739 events


## MAKE MODELS ----
# adapt
p1 <- ggplot(
  final,
  aes(x = adapt, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('adapt vs. microglianess') +
  ylab('microglianess') + xlab ('adapt') +
  theme_bw()

microglianess <- final$Micro.PVM
adapt <- final$adapt
m1 <- lm(microglianess ~ adapt)
summary(m1)

# average rate
p2 <- ggplot(
  final,
  aes(x = avg_rate, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('average rate vs. microglianess') +
  ylab('microglianess') + xlab ('average rate') +
  theme_bw()

avg_rate <- final$avg_rate
m2 <- lm(microglianess ~ avg_rate)
summary(m2)

# first isi
p3 <- ggplot(
  final,
  aes(x = first_isi, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('first isi vs. microglianess') +
  ylab('microglianess') + xlab ('first isi') +
  theme_bw()

first_isi <- final$first_isi
m3 <- lm(microglianess ~ first_isi)
summary(m3)

# isi cv
p4 <- ggplot(
  final,
  aes(x = isi_cv, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('isi cv vs. microglianess') +
  ylab('microglianess') + xlab ('isi cv') +
  theme_bw()

isi_cv <- final$isi_cv
m4 <- lm(microglianess ~ isi_cv)
summary(m4)

# latency
p5 <- ggplot(
  final,
  aes(x = latency, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('latency vs. microglianess') +
  ylab('microglianess') + xlab ('latency') +
  theme_bw()

latency <- final$latency
m5 <- lm(microglianess ~ latency)
summary(m5)

# mean isi
p6 <- ggplot(
  final,
  aes(x = mean_isi, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('mean isi vs. microglianess') +
  ylab('microglianess') + xlab ('mean isi') +
  theme_bw()

mean_isi <- final$mean_isi
m6 <- lm(microglianess ~ mean_isi)
summary(m6)

# median isi
p7 <- ggplot(
  final,
  aes(x = median_isi, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('median isi vs. microglianess') +
  ylab('microglianess') + xlab ('median isi') +
  theme_bw()

median_isi <- final$median_isi
m7 <- lm(microglianess ~ median_isi)
summary(m7)

# stim_amp 
p8 <- ggplot(
  final,
  aes(x = stim_amp, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('stim_amp vs. microglianess') +
  ylab('microglianess') + xlab ('stim_amp') +
  theme_bw()

stim_amp <- final$stim_amp
m8 <- lm(microglianess ~ stim_amp)
summary(m8)

# threshold_v 
p9 <- ggplot(
  final,
  aes(x = threshold_v, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('threshold_v vs. microglianess') +
  ylab('microglianess') + xlab ('threshold_v') +
  theme_bw()

threshold_v <- final$threshold_v
m9 <- lm(microglianess ~ threshold_v)
summary(m9)

# peak_v 
p10 <- ggplot(
  final,
  aes(x = peak_v, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('peak_v vs. microglianess') +
  ylab('microglianess') + xlab ('peak_v') +
  theme_bw()

peak_v <- final$peak_v
m10 <- lm(microglianess ~ peak_v)
summary(m10)

# trough_v 
p11 <- ggplot(
  final,
  aes(x = trough_v, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('trough_v vs. microglianess') +
  ylab('microglianess') + xlab ('trough_v') +
  theme_bw()

trough_v <- final$trough_v
m11 <- lm(microglianess ~ trough_v)
summary(m11)

# fast_trough_v
p12 <- ggplot(
  final,
  aes(x = fast_trough_v, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('fast_trough_v vs. microglianess') +
  ylab('microglianess') + xlab ('fast_trough_v') +
  theme_bw()

fast_trough_v <- final$fast_trough_v
m12 <- lm(microglianess ~ fast_trough_v)
summary(m12)

# adp_v 
p13 <- ggplot(
  final,
  aes(x = adp_v, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('adp_v vs. microglianess') +
  ylab('microglianess') + xlab ('adp_v') +
  theme_bw()

adp_v <- final$adp_v
m13 <- lm(microglianess ~ adp_v)
summary(m13)

# width 
p14 <- ggplot(
  final,
  aes(x = width, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('width vs. microglianess') +
  ylab('microglianess') + xlab ('width') +
  theme_bw()

width <- final$width
m14 <- lm(microglianess ~ width)
summary(m14)

# upstroke_downstroke_ratio 
p15 <- ggplot(
  final,
  aes(x = upstroke_downstroke_ratio, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('upstroke_downstroke_ratio vs. microglianess') +
  ylab('microglianess') + xlab ('upstroke_downstroke_ratio') +
  theme_bw()

upstroke_downstroke_ratio <- final$upstroke_downstroke_ratio
m15 <- lm(microglianess ~ upstroke_downstroke_ratio)
summary(m15)

# peak_t 
p16 <- ggplot(
  final,
  aes(x = peak_t, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('peak_t vs. microglianess') +
  ylab('microglianess') + xlab ('peak_t') +
  theme_bw()

peak_t <- final$peak_t
m16 <- lm(microglianess ~ peak_t)
summary(m16)

# fast_trough_t 
p17 <- ggplot(
  final,
  aes(x = fast_trough_t, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('fast_trough_t vs. microglianess') +
  ylab('microglianess') + xlab ('fast_trough_t') +
  theme_bw()

fast_trough_t <- final$fast_trough_t
m17 <- lm(microglianess ~ fast_trough_t)
summary(m17)

# trough_t 
p18 <- ggplot(
  final,
  aes(x = trough_t, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('trough_t vs. microglianess') +
  ylab('microglianess') + xlab ('trough_t') +
  theme_bw()

trough_t <- final$trough_t
m18 <- lm(microglianess ~ trough_t)
summary(m18)

# slow_trough_t 
p19 <- ggplot(
  final,
  aes(x = slow_trough_t, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('slow_trough_t vs. microglianess') +
  ylab('microglianess') + xlab ('slow_trough_t') +
  theme_bw()

slow_trough_t <- final$slow_trough_t
m19 <- lm(microglianess ~ slow_trough_t)
summary(m19)

# rheo_first_isi 
p20 <- ggplot(
  final,
  aes(x = rheo_first_isi, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('rheo_first_isi vs. microglianess') +
  ylab('microglianess') + xlab ('rheo_first_isi') +
  theme_bw()

rheo_first_isi <- final$rheo_first_isi
m20 <- lm(microglianess ~ rheo_first_isi)
summary(m20)

# v_baseline 
p21 <- ggplot(
  final,
  aes(x = v_baseline, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('v_baseline vs. microglianess') +
  ylab('microglianess') + xlab ('v_baseline') +
  theme_bw()

v_baseline <- final$v_baseline
m21 <- lm(microglianess ~ v_baseline)
summary(m21)

# rheobase_i 
p22 <- ggplot(
  final,
  aes(x = rheobase_i, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('rheobase_i vs. microglianess') +
  ylab('microglianess') + xlab ('rheobase_i') +
  theme_bw()

rheobase_i <- final$rheobase_i
m22 <- lm(microglianess ~ rheobase_i)
summary(m22)

# fi_fit_slope 
p23 <- ggplot(
  final,
  aes(x = fi_fit_slope, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('fi_fit_slope vs. microglianess') +
  ylab('microglianess') + xlab ('fi_fit_slope') +
  theme_bw()

fi_fit_slope <- final$fi_fit_slope
m23 <- lm(microglianess ~ fi_fit_slope)
summary(m23)

# sag 
p24 <- ggplot(
  final,
  aes(x = sag, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('sag vs. microglianess') +
  ylab('microglianess') + xlab ('sag') +
  theme_bw()

sag <- final$sag
m24 <- lm(microglianess ~ sag)
summary(m24)

# vm_for_sag 
p25 <- ggplot(
  final,
  aes(x = vm_for_sag, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('vm_for_sag vs. microglianess') +
  ylab('microglianess') + xlab ('vm_for_sag') +
  theme_bw()

vm_for_sag <- final$vm_for_sag
m25 <- lm(microglianess ~ vm_for_sag)
summary(m25)

# input_resistance 
p26 <- ggplot(
  final,
  aes(x = input_resistance, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('input_resistance vs. microglianess') +
  ylab('microglianess') + xlab ('input_resistance') +
  theme_bw()

input_resistance <- final$input_resistance
m26 <- lm(microglianess ~ input_resistance)
summary(m26)

# tau 
p27 <- ggplot(
  final,
  aes(x = tau, y = Micro.PVM)) +
  geom_smooth(method = "lm", se = F) + 
  geom_point() +
  ggtitle('tau vs. microglianess') +
  ylab('microglianess') + xlab ('tau') +
  theme_bw()

tau <- final$tau
m27 <- lm(microglianess ~ tau)
summary(m27)
