source("src/a-setup.R")

# get colony data
qdf.raw <- readxl::read_xlsx(path = file.path(data_dir, "colony-data_2024-08-08.xlsx"))
names(qdf.raw) <- make.names(names(qdf.raw))

qdf <- qdf.raw %>% 
  mutate(n.queens = Closed.Queen.Cells..Emergency. + Closed.Queen.Cells..Regular. + 
         Open.Queen.Cells..Emergency. +Open.Queen.Cells..Regular.,
         r.effort = n.queens/Workers,
         doy = yday(Date.found),
         yr = year(Date.found),
         yr.index = as.numeric(as.factor(yr)),
         doy.cont = doy + c(0, 365)[yr.index],
         mth = month(Date.found)) %>%
  select(date = Date.found, doy, doy.cont, mth, n.queens, r.effort)

ggplot(data = qdf, mapping = aes(x = date, y = r.effort)) +
  geom_point()

ggplot(data = qdf, mapping = aes(x = date, y = n.queens)) +
  geom_point()
 
month.summ <- mutate(qdf, mth = factor(mth, levels = 1:12)) %>%
  group_by(mth, .drop = FALSE) %>%
  summarise(mean.nq = mean(n.queens),
            mean.re = mean(r.effort))

ggplot(data = month.summ, mapping = aes(x = mth, y = mean.nq)) +
  geom_point()

ggplot(data = month.summ, mapping = aes(x = mth, y = mean.re)) +
  geom_point()
