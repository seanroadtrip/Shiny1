# ============================================================
# Restoration siting tool — South Island, New Zealand
#
# Draw a boundary, set operational constraints, and the app
# returns: susceptible land cover within the boundary, the
# portion restorable under those constraints, which indigenous
# species are environmentally suitable, and projected time to
# canopy closure for species combinations with growth data.
#
# All modelling is pre-computed. The app only clips and
# summarises rasters, so it stays fast and light.
#
# Save as app.R in the project root and run with shiny::runApp()
# ============================================================

library(shiny)
library(leaflet)
library(leaflet.extras)
library(terra)
library(sf)
library(DT)

# ------------------------------------------------------------
# Load pre-computed layers (lazy: values read only when cropped)
# ------------------------------------------------------------
setwd("C:/Users/LENOVO P15V G3/Desktop/R_Projects/SeanShiny")
surf_dir <- "suitability_surfaces"

susc_class  <- rast("susceptible_classes_100m.tif")
consensus   <- rast("treeline_consensus.tif")
dist_access <- rast("dist_to_access_100m.tif")
slope       <- rast("slope_100m.tif")

cell_km2 <- prod(res(susc_class)) / 1e6
cell_ha  <- prod(res(susc_class)) / 10000

lcdb_names <- c("15" = "Alpine Grass/Herbfield",  "41" = "Low Producing Grassland",
                "44" = "Depleted Grassland",      "51" = "Gorse and/or Broom",
                "55" = "Sub Alpine Shrubland",    "56" = "Mixed Exotic Shrubland",
                "58" = "Matagouri or Grey Scrub", "64" = "Forest - Harvested")

fits <- readRDS("cv_results.rds")
sp_names <- names(fits)
safe_name <- function(x) gsub(" ", "_", x)

# Model performance travels with each species so the interface
# can show how much confidence the envelope deserves
sp_perf <- do.call(rbind, lapply(fits, function(r) data.frame(
  species = r$species, auc = r$auc, boyce = r$boyce, omission = r$om_S05,
  n = r$n, stringsAsFactors = FALSE)))

conifers <- c("Pinus radiata", "Pinus contorta", "Pseudotsuga menziesii")

closure <- readRDS("closure_model_5sp.rds")

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Restoration siting tool — South Island"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      helpText("Draw a boundary on the map using the polygon or rectangle tool,",
               "then set the operational constraints below."),
      hr(),
      sliderInput("dist_km", "Maximum distance from road or track (km)",
                  min = 0.5, max = 8, value = 4, step = 0.5),
      sliderInput("slope_deg", "Maximum working slope (degrees)",
                  min = 5, max = 30, value = 14.4, step = 0.5),
      helpText(tags$small(
        "23 degrees (about 42 % grade) was used as a health and safety",
        "limit for volunteer planting crews on unimproved ground.")),
      hr(),
      sliderInput("suit_thresh", "Minimum environmental similarity for a species",
                  min = 0.01, max = 0.5, value = 0.05, step = 0.01),
      helpText(tags$small(
        "Lower values are more permissive. At 0.05, envelopes excluded",
        "about 21 % of sites where species were actually recorded,",
        "so a lower threshold may suit an advisory filter.")),
      hr(),
      actionButton("go", "Analyse", class = "btn-primary btn-block")
    ),
    
    mainPanel(
      width = 9,
      leafletOutput("map", height = 480),
      br(),
      tabsetPanel(
        tabPanel("Summary",   br(), uiOutput("summary_box"),
                 br(), DTOutput("cover_tbl")),
        tabPanel("Species",   br(), helpText(
          "Species whose environmental envelope covers the site, ranked by mean",
          "similarity. AUC and omission indicate how much confidence the envelope",
          "deserves: higher AUC means better discrimination, lower omission means",
          "the envelope less often excludes sites where the species does occur."),
          DTOutput("species_tbl")),
        tabPanel("Canopy closure", br(), uiOutput("closure_note"),
                 DTOutput("closure_tbl"))
      )
    )
  )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

server <- function(input, output, session) {
  
  output$map <- renderLeaflet({
    leaflet() |>
      addTiles(group = "Map") |>                                  # OpenStreetMap
      addProviderTiles("Esri.WorldImagery", group = "Satellite") |>
      addLayersControl(baseGroups = c("Map", "Satellite")) |>
      setView(lng = 171.0, lat = -43.5, zoom = 6) |>
      addDrawToolbar(
        targetGroup = "drawn", polylineOptions = FALSE, circleOptions = FALSE,
        markerOptions = FALSE, circleMarkerOptions = FALSE,
        editOptions = editToolbarOptions())
  })
  
  poly <- reactiveVal(NULL)
  
  observeEvent(input$map_draw_new_feature, {
    f <- input$map_draw_new_feature
    cds <- f$geometry$coordinates[[1]]
    m <- do.call(rbind, lapply(cds, function(p) c(p[[1]], p[[2]])))
    p <- st_sfc(st_polygon(list(m)), crs = 4326) |> st_transform(2193)
    poly(vect(st_sf(geometry = p)))
  })
  
  result <- eventReactive(input$go, {
    req(poly())
    p <- poly()
    
    d_m <- input$dist_km * 1000
    s_d <- input$slope_deg
    
    # Clip the constraint stack to the drawn boundary
    sc <- mask(crop(susc_class,  p), p)
    tl <- mask(crop(consensus,   p), p) >= 1
    da <- mask(crop(dist_access, p), p)
    sl <- mask(crop(slope,       p), p)
    
    n_susc <- global(!is.na(sc), "sum", na.rm = TRUE)[[1]]
    if (is.na(n_susc) || n_susc == 0)
      return(list(empty = TRUE, area_ha = expanse(p, unit = "ha")[1]))
    
    rest <- !is.na(sc) & !tl & (da <= d_m) & (sl <= s_d)
    
    cover <- freq(mask(sc, rest, maskvalues = c(FALSE, NA)))
    cover_all <- freq(sc)
    
    tab <- data.frame(
      Class = lcdb_names[as.character(cover_all$value)],
      Susceptible_ha = round(cover_all$count * cell_ha),
      Restorable_ha  = 0, stringsAsFactors = FALSE)
    if (nrow(cover))
      tab$Restorable_ha[match(cover$value, cover_all$value)] <-
      round(cover$count * cell_ha)
    tab$Percent <- round(100 * tab$Restorable_ha / tab$Susceptible_ha, 1)
    tab <- tab[order(-tab$Susceptible_ha), ]
    
    # Mean environmental similarity for each species within the boundary
    suit <- do.call(rbind, lapply(sp_names, function(sp) {
      r <- rast(file.path(surf_dir, paste0(safe_name(sp), ".tif")))
      v <- mask(crop(r, p), p)
      data.frame(species = sp,
                 mean_suit = global(v, "mean", na.rm = TRUE)[[1]],
                 pct_above = 100 * global(v >= input$suit_thresh, "mean",
                                          na.rm = TRUE)[[1]])
    }))
    suit <- merge(suit, sp_perf, by = "species")
    suit <- suit[order(-suit$pct_above), ]
    
    list(empty = FALSE,
         area_ha    = expanse(p, unit = "ha")[1],
         susc_ha    = round(n_susc * cell_ha),
         rest_ha    = round(global(rest, "sum", na.rm = TRUE)[[1]] * cell_ha),
         above_tl   = round(global(!is.na(sc) & tl, "sum", na.rm = TRUE)[[1]] * cell_ha),
         too_steep  = round(global(!is.na(sc) & !tl & sl > s_d, "sum", na.rm = TRUE)[[1]] * cell_ha),
         too_far    = round(global(!is.na(sc) & !tl & sl <= s_d & da > d_m, "sum", na.rm = TRUE)[[1]] * cell_ha),
         cover = tab, suit = suit)
  })
  
  output$summary_box <- renderUI({
    r <- result()
    if (r$empty)
      return(div(class = "alert alert-warning",
                 sprintf("No invasion-susceptible land cover within this %s ha boundary.",
                         format(round(r$area_ha), big.mark = ","))))
    f <- function(x) format(x, big.mark = ",")
    tagList(
      h4(sprintf("%s ha boundary", f(round(r$area_ha)))),
      tags$ul(
        tags$li(sprintf("Invasion-susceptible land cover: %s ha", f(r$susc_ha))),
        tags$li(tags$b(sprintf("Restorable under these constraints: %s ha (%.0f %% of susceptible)",
                               f(r$rest_ha), 100 * r$rest_ha / r$susc_ha))),
        tags$li(sprintf("Excluded — above treeline: %s ha", f(r$above_tl))),
        tags$li(sprintf("Excluded — too steep: %s ha", f(r$too_steep))),
        tags$li(sprintf("Excluded — too far from access: %s ha", f(r$too_far)))))
  })
  
  output$cover_tbl <- renderDT({
    r <- result(); req(!r$empty)
    datatable(r$cover, rownames = FALSE, options = list(dom = "t", paging = FALSE),
              colnames = c("Land cover class", "Susceptible (ha)",
                           "Restorable (ha)", "% restorable"))
  })
  
  output$species_tbl <- renderDT({
    r <- result(); req(!r$empty)
    d <- r$suit
    d$Type <- ifelse(d$species %in% conifers, "Wilding conifer", "Indigenous")
    out <- data.frame(
      Species = d$species, Type = d$Type,
      `Site coverage (%)` = round(d$pct_above, 1),
      `Mean similarity`   = round(d$mean_suit, 3),
      AUC = round(d$auc, 2), `Omission (%)` = round(d$omission, 1),
      Records = d$n, check.names = FALSE)
    datatable(out, rownames = FALSE, options = list(pageLength = 18, dom = "t")) |>
      formatStyle("Type", target = "row",
                  backgroundColor = styleEqual("Wilding conifer", "#f2dede"))
  })
  
  output$closure_note <- renderUI({
    r <- result(); req(!r$empty)
    ok <- intersect(closure$species |> (\(s) closure$full_name[s])(),
                    r$suit$species[r$suit$pct_above >= 50])
    tagList(
      helpText(
        "Projected years to reach LAI 2.5, the canopy closure threshold associated with",
        "resistance to wilding conifer establishment. Estimates derive from 185 plots at",
        sprintf("43 restoration sites aged %.1f to %.0f years, and are indicative rather than predictive.",
                closure$age_range[["min"]], closure$age_range[["max"]]),
        "Only combinations observed in at least five plots are shown. Species are treated",
        "as present or absent; planting density is not accounted for."),
      if (!length(ok)) div(class = "alert alert-info",
                           "None of the five species with growth data is environmentally suitable across at least half of this site.")
    )
  })
  
  output$closure_tbl <- renderDT({
    r <- result(); req(!r$empty)
    suitable <- r$suit$species[r$suit$pct_above >= 50]
    lk <- closure$lookup[closure$lookup$n_plots >= 5, ]
    
    # Keep only combinations whose species are all suitable at this site
    keep <- sapply(strsplit(lk$species, "\\+"), function(codes) {
      if (identical(codes, "none")) return(TRUE)
      all(closure$full_name[codes] %in% suitable)
    })
    lk <- lk[keep, c("assemblage", "age_closure", "n_plots", "extrapolated")]
    
    datatable(lk, rownames = FALSE,
              options = list(pageLength = 15, order = list(list(1, "asc"))),
              colnames = c("Species assemblage", "Years to LAI 2.5",
                           "Supporting plots", "Beyond observed ages"))
  })
}

shinyApp(ui, server)


rsconnect::writeManifest()



