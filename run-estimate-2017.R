# Estimate daily and total herring run size using DMF methodology
#
# Uses stratified 2-way random sampling (St2WRS)
# with three 4-hour periods (0700-1100, 1100-1500, 1500-1900)
#
# Jeffrey D Walker, PhD <jeff@walkerenvres.com>
# Walker Environmental Research, LLC <https://walkerenvres.com>
# March 26, 2018
#
# Ref:
#   Nelson, 2006
#   A Guide to Statistical Sampling for the Estimation of River Herring Run Size Using Visual Counts
#   MA DMF Technical Report TR-25
#   http://www.mass.gov/eea/docs/dfg/dmf/publications/tr-25.pdf

library(tidyverse)
library(lubridate)


# setup -------------------------------------------------------------------

# load volunteer count data
df_counts_all <- read_csv(
  "volunteer-counts-2017.csv",
  col_types = cols(
    Date = col_character(),
    Collector = col_character(),
    Start_Time = col_time(format = ""),
    End_Time = col_time(format = ""),
    Count = col_integer()
  )
) %>%
  mutate(
    start = mdy_hms(paste(Date, Start_Time, sep = " ")),
    end = mdy_hms(paste(Date, End_Time, sep = " ")),
    date = mdy(Date)
  ) %>%
  select(
    collector = Collector,
    date,
    start,
    end,
    count = Count
  )

# list of volunteer counts to exclude (end timestamps)
exclude_ends <- ymd_hms(
  "2017-05-27 07:51:00", # ?? not sure why this one was excluded
  "2017-05-28 10:50:00", # overlaps other count
  "2017-06-06 18:51:00", # overlaps other count
  "2017-06-24 07:00:00"  # too early
)

# filter counts (exclude too early, too late), and assign periods
df_counts <- df_counts_all %>%
  filter(
    month(date) < 7,
    hour(start) < 19,
    hour(end) >= 7,
    !end %in% exclude_ends,
    !(hour(start) != hour(end) & hour(end) %in% c(7, 11, 15, 19) & minute(end) > 0) # exclude any 10-min periods that cross from one 4-hour period to the next
  ) %>%
  mutate(
    hour = hour(start),
    period = case_when(
      hour >= 7 & hour < 11 ~ 1,
      hour >= 11 & hour < 15 ~ 2,
      hour >= 15 & hour < 19 ~ 3,
      TRUE ~ NA_real_
    )
  ) %>%
  rename(
    y_kpi = count
  ) %>%
  select(date, start, end, hour, period, y_kpi)

# daily summary statistics, note that these are NOT used for total run estimates
df_daily_stats <- df_counts %>%
  group_by(date) %>%
  summarise(
    mean = mean(y_kpi),
    sd = sd(y_kpi),
    n = n(),
    se = sd / sqrt(n)
  )


# calculations ------------------------------------------------------------

df_period <- df_counts %>%
  group_by(date, period) %>%
  summarise(
    y_kp = mean(y_kpi),            # mean count per period
    n_kp = n(),
    s_kp = coalesce(sd(y_kpi), 0), # st. dev of period count (Eq 13)
    N_kp = 6 * 4
  ) %>%
  ungroup() %>%
  mutate(
    Y_kp = N_kp * y_kp,                              # estimated total coun per period (Eq 12)
    varY_kp = N_kp * (N_kp - n_kp) * s_kp^2 / n_kp,  # variance of total count per period (Eq 13)
    a_kp = N_kp * (N_kp - n_kp) / n_kp               # coefficient used to estimate deg. of freedom (Eq 9)
  )

df_day <- df_period %>%
  group_by(date) %>%
  summarise(
    n_p = n(),                                      # Number of periods per day
    n_k = sum(n_kp),                                # Number of counts per day
    N_k = sum(N_kp),                                # Number of 10-min intervals per day (NOTE! only includes daily hours, 0700-1900)
    Y_k = sum(Y_kp),                                # Daily estimated run (Eq 12)
    varY_k = sum(varY_kp),                          # Variance of daily estimated run (Eq 13)
    se_k = sqrt(varY_k),
    df_num = sum(a_kp * s_kp ^ 2),                  # numerator of Eq 9
    df_den = sum((a_kp * s_kp ^ 2) ^ 2 / (n_kp - 1), na.rm = TRUE), # denominator of Eq 9
    df = df_num ^ 2 / df_den,                       # degrees of freedom (Eq 9)
    df = round(df),
    t_star = map_dbl(df, ~ qt(0.975, df = .)),      # t distribution statistic, two-sided 95% CI
    df = coalesce(df, 0),
    t_star = coalesce(t_star, 0),
    ci_lower = Y_k - t_star * sqrt(varY_k),         # Lower 95% CI Bound (Eq 14)
    ci_upper = Y_k + t_star * sqrt(varY_k)          # Upper 95% CI Bound (Eq 14)
  ) %>%
  mutate(
    Y_k = round(Y_k),
    ci_lower = round(ci_lower),
    ci_upper = round(ci_upper)
  )

df_run <- df_day %>%
  summarise(
    N = sum(N_k),
    total = sum(Y_k),                            # Total run estimate
    se = sqrt(sum(varY_k)),                      # Std. error of total run estimate
    df = sum(df_num) ^ 2 / sum(df_den),          # Degrees of freedom (Eq 9)
    df = round(df),
    t_star = map_dbl(df, ~ qt(0.975, df = .)),   # t statistic
    ci_lower = total - t_star * se,              # Lower 95% CI Bound (Eq 14)
    ci_upper = total + t_star * se               # Upper 95% CI Bound (Eq 14)
  )


# results -----------------------------------------------------------------

print(df_counts) # individual counts
print(df_period) # period statistics
print(df_day)    # daily statistics
print(df_run)    # total run statistics

# plot daily run estimates with 95% CI
df_day %>%
  ggplot(aes(date, Y_k)) +
  geom_col() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper)) +
  ylim(0, NA) +
  labs(
    x = "Date",
    y = "Estimate Daily Run (95% CI)"
  )


# export ------------------------------------------------------------------

df_day %>%
  write_csv("run-estimate-2017-daily.csv", na = "")
df_run %>%
  write_csv("run-estimate-2017-total.csv", na = "")
