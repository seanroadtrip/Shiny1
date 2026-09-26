# ============================================================
# Restoration siting tool — South Island, New Zealand
#
# Draw a boundary, set operational constraints, and the app
# returns: susceptible land cover, the portion restorable under
# those constraints, which species are environmentally suitable,
# and projected time to canopy closure for species combinations
# with growth data.
#
# All modelling is pre-computed; the app clips and summarises.
#
# Folder layout:
#   app.R
#   GIS_files/suitability_surfaces/  (one .tif per species)
#   GIS_files/susceptible_classes_100m.tif
#   GIS_files/treeline_consensus.tif
#   GIS_files/dist_to_access_100m.tif
#   GIS_files/slope_100m.tif
#   outputs/cv_results.rds
#   outputs/closure_model_5sp.rds
# ============================================================

surf_dir <- "suitability_surfaces"

susc_class  <- rast("susceptible_classes_100m.tif")
consensus   <- rast("treeline_consensus.tif")
dist_access <- rast("dist_to_access_100m.tif")
slope       <- rast("slope_100m.tif")

cell_ha <- prod(res(susc_class)) / 10000

lcdb_names <- c("15" = "Alpine Grass/Herbfield",  "41" = "Low Producing Grassland",
                "44" = "Depleted Grassland",      "51" = "Gorse and/or Broom",
                "55" = "Sub Alpine Shrubland",    "56" = "Mixed Exotic Shrubland",
                "58" = "Matagouri or Grey Scrub", "64" = "Forest - Harvested")
lcdb_cols  <- c("41" = "#A3D400", "55" = "#B8AB6A", "15" = "#ABCD66",
                "44" = "#D2D25A", "51" = "#7D690F", "58" = "#D4CDAE",
                "64" = "#A1AD61", "56" = "#C4BB89")

fits <- readRDS("cv_results.rds")
sp_names  <- names(fits)
safe_name <- function(x) gsub(" ", "_", x)

sp_perf <- do.call(rbind, lapply(fits, function(r) data.frame(
  species = r$species, auc = r$auc, omission = r$om_S05, n = r$n,
  stringsAsFactors = FALSE)))

conifers <- c("Pinus radiata", "Pinus contorta", "Pseudotsuga menziesii")

closure <- readRDS("closure_model_5sp.rds")

# Species with individually significant coefficients in the best
# supported model; the other three appeared in competitive models
# but were not significant on their own.
sig_species <- c("COPROB", "KUNERI")

support_for <- function(code_string) {
  if (identical(code_string, "none")) return(list(w = NA_real_, level = "Baseline"))
  codes <- strsplit(code_string, "\\+")[[1]]
  w <- min(closure$akaike_weight[codes])
  lvl <- if (all(codes %in% sig_species)) "Strong"
  else if (w >= 0.5) "Moderate" else "Weak"
  list(w = w, level = lvl)
}

# ------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Restoration siting tool — South Island"),
  
  sidebarLayout(
    sidebarPanel(width = 3,
                 strong("1. Draw your site"),
                 helpText("Use the polygon or rectangle tool on the map."),
                 hr(),
                 strong("2. Set your constraints"),
                 sliderInput("dist_km", "Furthest crews will carry plants (km)",
                             min = 0.5, max = 8, value = 4, step = 0.5),
                 sliderInput("slope_deg", "Steepest ground crews will work (degrees)",
                             min = 5, max = 30, value = 14.4, step = 0.5),
                 helpText(tags$small("23 degrees is roughly a 42 % grade, used here as a",
                                     "health and safety limit for volunteer crews on",
                                     "unimproved ground.")),
                 sliderInput("suit_thresh", "How strict should the species filter be?",
                             min = 0.01, max = 0.5, value = 0.05, step = 0.01),
                 helpText(tags$small("Lower values include more species. At 0.05 the",
                                     "envelopes excluded roughly one site in five where a",
                                     "species was in fact recorded, so being permissive is",
                                     "reasonable.")),
                 hr(),
                 actionButton("go", "Analyse", class = "btn-primary btn-block"),
                 hr(),
                 uiOutput("sp_picker")
    ),
    
    mainPanel(width = 9,
              leafletOutput("map", height = 480), br(),
              tabsetPanel(
                tabPanel("Summary", br(),
                         helpText("Land cover classes vulnerable to wilding conifer invasion within",
                                  "your site, and how much of that area can realistically be planted",
                                  "given the constraints you set. Shaded areas on the map show where."),
                         uiOutput("summary_box"), br(), DTOutput("cover_tbl")),
                
                tabPanel("Species", br(),
                         helpText(tags$b("Which species suit this site?"),
                                  "Each species is ranked by how much of your site falls inside the",
                                  "environmental conditions where it has been recorded across the South",
                                  "Island. Site coverage is the percentage of your site inside that",
                                  "envelope."),
                         helpText(tags$small(
                           tags$b("Reading the confidence columns: "),
                           "AUC measures how well the model distinguishes places the species grows",
                           "from places botanists have surveyed generally; 0.5 is no better than",
                           "chance and 1.0 is perfect. Omission is how often the model wrongly",
                           "excluded sites where the species was actually recorded, so lower is",
                           "better. Widespread species tend to have low AUC and low omission: poor",
                           "at discriminating, but reliable as a filter. Wilding conifers are shown",
                           "for reference and are shaded red.")),
                         DTOutput("species_tbl")),
                
                tabPanel("Canopy closure", br(), uiOutput("closure_note"),
                         DTOutput("closure_tbl"))
              )
    )
  )
)

# ------------------------------------------------------------
server <- function(input, output, session) {
  
  output$map <- renderLeaflet({
    leaflet() |>
      addTiles(group = "Map") |>
      addProviderTiles("Esri.WorldImagery", group = "Satellite") |>
      addLayersControl(baseGroups = c("Map", "Satellite")) |>
      setView(lng = 171.0, lat = -43.5, zoom = 6) |>
      addDrawToolbar(targetGroup = "drawn", polylineOptions = FALSE,
                     circleOptions = FALSE, markerOptions = FALSE,
                     circleMarkerOptions = FALSE,
                     editOptions = editToolbarOptions())
  })
  
  poly <- reactiveVal(NULL)
  
  observeEvent(input$map_draw_new_feature, {
    cds <- input$map_draw_new_feature$geometry$coordinates[[1]]
    m <- do.call(rbind, lapply(cds, function(p) c(p[[1]], p[[2]])))
    poly(vect(st_sf(geometry = st_transform(
      st_sfc(st_polygon(list(m)), crs = 4326), 2193))))
  })
  
  result <- eventReactive(input$go, {
    req(poly()); p <- poly()
    d_m <- input$dist_km * 1000; s_d <- input$slope_deg
    
    sc <- mask(crop(susc_class,  p), p)
    tl <- mask(crop(consensus,   p), p) >= 1
    da <- mask(crop(dist_access, p), p)
    sl <- mask(crop(slope,       p), p)
    
    n_susc <- global(!is.na(sc), "sum", na.rm = TRUE)[[1]]
    if (is.na(n_susc) || n_susc == 0)
      return(list(empty = TRUE, area_ha = expanse(p, unit = "ha")[1]))
    
    rest <- !is.na(sc) & !tl & (da <= d_m) & (sl <= s_d)
    
    rest_cls  <- mask(sc, rest, maskvalues = c(FALSE, NA))
    rest_poly <- if (expanse(p, unit = "ha")[1] < 500000) {
      tryCatch(st_transform(st_make_valid(st_as_sf(
        as.polygons(rest_cls, dissolve = TRUE, na.rm = TRUE))), 4326),
        error = function(e) NULL)
    } else NULL
    
    cover_all <- freq(sc); cover <- freq(rest_cls)
    tab <- data.frame(Class = lcdb_names[as.character(cover_all$value)],
                      Susceptible_ha = round(cover_all$count * cell_ha),
                      Restorable_ha = 0, stringsAsFactors = FALSE)
    if (nrow(cover))
      tab$Restorable_ha[match(cover$value, cover_all$value)] <-
      round(cover$count * cell_ha)
    tab$Percent <- round(100 * tab$Restorable_ha / tab$Susceptible_ha, 1)
    tab <- tab[order(-tab$Susceptible_ha), ]
    
    suit <- withProgress(message = "Assessing species", value = 0, {
      do.call(rbind, lapply(seq_along(sp_names), function(i) {
        sp <- sp_names[i]
        incProgress(1 / length(sp_names), detail = sp)
        f <- file.path(surf_dir, paste0(safe_name(sp), ".tif"))
        if (!file.exists(f)) return(NULL)
        v <- mask(crop(rast(f), p), p)
        d <- data.frame(species = sp,
                        pct_above = 100 * global(v >= input$suit_thresh, "mean",
                                                 na.rm = TRUE)[[1]],
                        mean_suit = global(v, "mean", na.rm = TRUE)[[1]])
        rm(v); d
      }))
    })
    
    validate(need(!is.null(suit) && nrow(suit) > 0, paste(
      "No suitability surfaces found in", surf_dir,
      "- check the .tif files are in the app folder.")))
    
    suit <- merge(suit, sp_perf, by = "species")
    suit <- suit[order(-suit$pct_above), ]
    
    list(empty = FALSE, area_ha = expanse(p, unit = "ha")[1],
         susc_ha = round(n_susc * cell_ha),
         rest_ha = round(global(rest, "sum", na.rm = TRUE)[[1]] * cell_ha),
         above_tl  = round(global(!is.na(sc) & tl, "sum", na.rm = TRUE)[[1]] * cell_ha),
         too_steep = round(global(!is.na(sc) & !tl & sl > s_d, "sum", na.rm = TRUE)[[1]] * cell_ha),
         too_far   = round(global(!is.na(sc) & !tl & sl <= s_d & da > d_m, "sum", na.rm = TRUE)[[1]] * cell_ha),
         cover = tab, suit = suit, rest_poly = rest_poly)
  })
  
  # Restorable area always drawn
  observeEvent(result(), {
    r <- result()
    proxy <- leafletProxy("map") |> clearGroup("Restorable")
    if (!r$empty && !is.null(r$rest_poly) && nrow(r$rest_poly)) {
      cls <- as.character(r$rest_poly[[1]])
      proxy |> addPolygons(data = r$rest_poly, group = "Restorable",
                           fillColor = unname(lcdb_cols[cls]), fillOpacity = 0.65,
                           color = "#1B4332", weight = 0.6,
                           label = unname(lcdb_names[cls]))
    }
  })
  
  # ---- Species picker for the closure tab ----
  output$sp_picker <- renderUI({
    r <- result(); req(!r$empty)
    codes <- closure$species
    nm    <- unname(closure$full_name[codes])
    cov   <- round(r$suit$pct_above[match(nm, r$suit$species)])
    tagList(
      strong("3. Choose your planting mix"),
      helpText(tags$small(
        "Growth data exists for these five species only. Percentages show how much",
        "of your site suits each. Tick the ones you could realistically source and",
        "plant; the Canopy closure tab will show how long each combination takes.")),
      checkboxGroupInput("plant_sp", NULL,
                         choiceNames  = sprintf("%s (%d %% of site)", nm, cov),
                         choiceValues = codes,
                         selected     = codes[!is.na(cov) & cov >= 25]))
  })
  
  output$summary_box <- renderUI({
    r <- result()
    if (r$empty) return(div(class = "alert alert-warning",
                            sprintf("No invasion-susceptible land cover within this %s ha boundary.",
                                    format(round(r$area_ha), big.mark = ","))))
    f <- function(x) format(x, big.mark = ",")
    tagList(h4(sprintf("%s ha site", f(round(r$area_ha)))),
            tags$ul(
              tags$li(sprintf("Vulnerable to wilding conifer invasion: %s ha", f(r$susc_ha))),
              tags$li(tags$b(sprintf("Plantable under your constraints: %s ha (%.0f %%)",
                                     f(r$rest_ha), 100 * r$rest_ha / r$susc_ha))),
              tags$li(sprintf("Ruled out, above the treeline: %s ha", f(r$above_tl))),
              tags$li(sprintf("Ruled out, too steep: %s ha", f(r$too_steep))),
              tags$li(sprintf("Ruled out, too far from access: %s ha", f(r$too_far)))),
            helpText(tags$small("Each area is counted once, against the first constraint",
                                "it failed: treeline, then slope, then distance.")))
  })
  
  output$cover_tbl <- renderDT({
    r <- result(); req(!r$empty)
    datatable(r$cover, rownames = FALSE, options = list(dom = "t", paging = FALSE),
              colnames = c("Land cover class", "Vulnerable (ha)",
                           "Plantable (ha)", "% plantable"))
  })
  
  output$species_tbl <- renderDT({
    r <- result(); req(!r$empty); d <- r$suit
    out <- data.frame(Species = d$species,
                      Type = ifelse(d$species %in% conifers, "Wilding conifer", "Indigenous"),
                      `Site coverage (%)` = round(d$pct_above, 1),
                      `Mean similarity` = round(d$mean_suit, 3),
                      AUC = round(d$auc, 2), `Omission (%)` = round(d$omission, 1),
                      `Records used` = d$n, check.names = FALSE)
    datatable(out, rownames = FALSE, options = list(pageLength = 18, dom = "t")) |>
      formatStyle("Type", target = "row",
                  backgroundColor = styleEqual("Wilding conifer", "#f2dede"))
  })
  
  output$closure_note <- renderUI({
    tagList(
      helpText(tags$b("How long until the canopy closes?"),
               "Once leaf area index reaches about 2.5, a planted stand shades the ground",
               "enough to resist wilding conifer seedlings establishing. The table shows how",
               "long each combination of the species you ticked is projected to take."),
      helpText(tags$small(
        sprintf("Projections come from 185 plots at 43 South Island restoration sites aged %.1f to %.0f years. ",
                closure$age_range[["min"]], closure$age_range[["max"]]),
        "They are indicative, not predictive: species are treated as present or absent,",
        "planting density is not accounted for, and combinations recorded in fewer than",
        "five plots are not shown.")),
      helpText(tags$small(tags$b("Support: "),
                          "Strong means every species in the combination had a statistically significant",
                          "effect on leaf area index. Moderate means every species appeared in most of the",
                          "competitive models but was not individually significant. Weak means at least one",
                          "species had limited support in the data. Treat Weak rows as rough guidance only.")))
  })
  
  output$closure_tbl <- renderDT({
    r <- result(); req(!r$empty)
    chosen <- input$plant_sp
    lk <- closure$lookup[closure$lookup$n_plots >= 5, ]
    
    keep <- sapply(strsplit(lk$species, "\\+"), function(codes)
      identical(codes, "none") || all(codes %in% chosen))
    lk <- lk[keep, ]
    
    if (!nrow(lk)) return(datatable(data.frame(
      Note = "Tick at least one species in the sidebar to see closure projections."),
      rownames = FALSE, options = list(dom = "t")))
    
    sup <- lapply(lk$species, support_for)
    out <- data.frame(Assemblage = lk$assemblage,
                      `Years to canopy closure` = lk$age_closure,
                      `Plots supporting this mix` = lk$n_plots,
                      Support = sapply(sup, `[[`, "level"),
                      check.names = FALSE)
    
    datatable(out, rownames = FALSE,
              options = list(pageLength = 15, order = list(list(1, "asc")))) |>
      formatStyle("Support", backgroundColor = styleEqual(
        c("Strong", "Moderate", "Weak"), c("#D8F3DC", "#FFF3CD", "#F8D7DA")))
  })
}

shinyApp(ui, server)
