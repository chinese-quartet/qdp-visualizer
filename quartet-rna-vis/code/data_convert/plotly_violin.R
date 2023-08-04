library(plotly)

dt <- readRDS('./rds/web/corfpkm.rds')
fig <- dt %>%
  plot_ly(type = 'violin')  %>%
  add_trace(
    x = ~batch_type[dt$libraryPrepA == 'PolyA'],
    y = ~Cor[dt$libraryPrepA == 'PolyA'],
    legendgroup = 'PolyA',
    scalegroup = 'PolyA',
    name = 'PolyA',
    side = 'negative',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = F
    ),
    color = I("blue")) %>%
  add_trace(
    x = ~batch_type[dt$libraryPrepA == 'RiboZero'],
    y = ~Cor[dt$libraryPrepA == 'RiboZero'],
    legendgroup = 'RiboZero',
    scalegroup = 'RiboZero',
    name = 'RiboZero',
    side = 'positive',
    box = list(
      visible = T
    ),
    meanline = list(
      visible = F
    ),
    color = I("orange")
  )  %>%
  layout(
    xaxis = list(
      title = ""  
    ),
    yaxis = list(
      title = "",
      zeroline = F
    ),
    violingap = 0,
    violingroupgap = 0,
    violinmode = 'overlay'
  )
  

fig

