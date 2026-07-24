 is a no-code analysis of profiler.Rprofvis.xml.

What the profile represents

┌─────────────────────┬──────────────────────────────────────────┐
│       Metric        │                  Value                   │
├─────────────────────┼──────────────────────────────────────────┤
│ Sampled ticks       │ 2,159                                    │
├─────────────────────┼──────────────────────────────────────────┤
│ Estimated wall time │ ~21.6 s (default profvis 10 ms interval) │
├─────────────────────┼──────────────────────────────────────────┤
│ Total stack frames  │ 349,197 (~162 frames per sample)         │
├─────────────────────┼──────────────────────────────────────────┤
│ flushReact cycles   │ 1,962                                    │
└─────────────────────┴──────────────────────────────────────────┘

Most of the runtime is not in data manipulation — it is in rendering, reactive invalidation, and repeatedly converting the same ggplots to plotly objects.

Hottest project-level functions (wall-clock tick coverage)

These are the CySA functions that appear most often in active sample frames:

┌──────────────────────┬────────────┬───────────────────────────────────────────┐
│       Function       │ % of ticks │               What it does                │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ quiet_ggplotly       │ ~22.5 %    │ Converts every ggplot to plotly           │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ somPlot              │ ~12.5 %    │ Builds + converts the 6 SOM scatter plots │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ VlnPlot2 / vlnPlot   │ ~11.8 %    │ Violin plots                              │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ umap::umap           │ ~3.7 %     │ UMAP dimension reduction                  │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ theme_minimal        │ ~2.5 %     │ Theme reconstruction                      │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ drawProjection       │ ~2.1 %     │ Projection overlays for SOM plots         │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ Rtsne::Rtsne         │ ~2.6 %     │ t-SNE computation                         │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ plotCytoScatter      │ ~1.1 %     │ Scatter density plots                     │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ upsetPlotFunc        │ ~0.7 %     │ ComplexHeatmap UpSet                      │
├──────────────────────┼────────────┼───────────────────────────────────────────┤
│ .buildFlowSOMPiePlot │ ~0.6 %     │ FlowSOM pie overlay                       │
└──────────────────────┴────────────┴───────────────────────────────────────────┘

Key rendering costs

The profile shows that rendering internals dominate:

- ggplotly / plotly_build / gg2list: appear together in ~22–25 % of ticks
- grid.draw / drawGTree: ~2.8 % of ticks, but heavily tied to the static outputs (VlnPlot2, flowSOMPie, somRaster, CountBar, PercentBar)
- ggplot_build, ggplot_gtable, calc_element, draw_panels, draw_panel: together ~6–7 % of ticks

So roughly two thirds of the profile is spent inside ggplot/plotly/grid rendering machinery, often re-doing work that has already been done.

Most important observations

1. The same plots are converted to plotly over and over
  - quiet_ggplotly is active in almost a quarter of all ticks.
  - The six somData outputs each run ~52–68 times and each call quiet_ggplotly.
  - output$umap, output$tsne, output$pca, output$somRasterSelect, and output$scatter also call it every render.
2. Outputs co-fire in the same flush cycles
  - output$umap and output$tsne fire together in 178 cycles.
  - The six somData plots fire in cascades.
  - output$VlnPlot2 and output$VlnPlot fire together in 85 cycles.
  - This means one user interaction can trigger 5–8 full plot rebuilds + plotly conversions.
3. Some outputs render even when nobody can see them
  - Several plot boxes are commented with req(input$violinBox)/req(input$upsetBox) but those guards are disabled.
  - VlnPlot2 runs 254 times, VlnPlot 141 times, UpSet 59 times.
4. Expensive dimension reductions are cached but still re-trigger
  - umap::umap and Rtsne::Rtsne are wrapped in bindCache, but they still appear in 3–4 % of ticks, suggesting the cache keys are not as stable as they could be or the reactives depend on too many transient inputs.
5. Static outputs still pay ggplot-to-plotly costs indirectly
  - VlnPlot2, flowSOMPie, somRaster, CountBar, and PercentBar use renderPlot, so they avoid ggplotly, but they still spend a lot of time in grid.draw, print.ggplot, and dev.off.
6. Theme objects are rebuilt every plot
  - theme_minimal() and theme_cowplot() appear in ~2.5 % and ~1 % of ticks respectively, even though the theme never changes.
7. The scatter plot adds 10,000 invisible trace points
  - output$scatter builds a 100×100 selection grid (xpoints, ypoints) and calls add_trace every render. That trace is not cheap.
8. drawProjection and somPlot redo the same joins
  - drawProjection re-orders df and re-attaches SOM_stats on every render.
  - somPlot re-runs quiet_ggplotly even when only the selected nodes changed.
9. Bar plots do unnecessary reshaping
  - countBarPlotFunc and PercentBarPlotFunc use reshape2::melt and dplyr::left_join/grouping on every render; the underlying clusterPatientTable and experiment_info rarely change.
10. A leftover <reactive:plotObj> consumes many samples
  - It appears in 994 frames and is a leaf frame, suggesting an old cached/observer reactive is still executing and doing nothing useful.

Proposed improvements (no code changes yet)

Grouped by expected impact.

High impact — reduce rendering work

1. Convert fewer plots to plotly
  - Keep the SOM mini-scatter plots (somData1–somData6) as static renderPlot outputs unless they genuinely need plotly selection. They are currently rebuilt and converted to plotly constantly.
  - Do the same for VlnPlot2, flowSOMPie, and UpSet if their interactivity is not required.
2. Cache the plotly object, not just the ggplot
  - The somBasePlots cache only the ggplot. Add a second layer that caches the result of quiet_ggplotly(...) keyed on the same inputs plus rs/colorbyGroups.
  - Similarly cache the plotly results for tsnePlot, umapPlot, pcaPlot, and somRasterSelect.
3. Combine the six SOM plots into one subplot
  - Instead of six independent renderPlotly outputs, build one plotly::subplot or faceted ggplot and convert once. This removes five ggplotly conversions and five sets of observers/event registrations per update.
4. Suspend outputs that are not visible
  - Re-enable req(input$violinBox) / req(input$upsetBox) guards, or use shinyjs::hidden + outputOptions(suspendWhenHidden = TRUE) so collapsed boxes do not re-render.
  - Apply the same to the FlowSOM stars/pie boxes when fsom is large.
5. Replace the invisible 10k-point selection trace
  - In output$scatter, use a plotly layout(shapes = list(...)) rectangle or the built-in dragmode = "select" instead of adding 10,000 transparent markers.

Medium impact — avoid repeated computation

6. Precompute and reuse theme objects
  - Create base_theme_minimal <- theme_minimal() + theme(...) and base_theme_cowplot <- theme_cowplot(...) + theme(...) once, then add them to plots instead of calling the constructors each time.
7. Stabilize dimension-reduction cache keys
  - Make sure tsne() and umap() reactives only depend on the columns, perplexity, and n_neighbors. Move selectedUpdate2() and triggerRedraw() out of the reactive body so they do not invalidate the cached result.
8. Precompute drawProjection data
  - Join SOM_stats into the projection data frame once when dimSelection changes, not on every render. Store row order.
9. Cache UpSet combinatorics
  - ComplexHeatmap::make_comb_mat and the two HeatmapAnnotation objects can be cached on upsetSelection and outputList so the UpSet plot only redraws when those actually change.
10. Memoize bar-plot data
  - Precompute the melted dfm for count/percentage bars whenever rs, cst, or groupsInput change, rather than building it inside the render function.

Cleanup

11. Remove the stale <reactive:plotObj>
  - Find and delete or rename the leftover plotObj reactive that is executing ~1,000 times.
12. Profile a focused interaction
  - The current profile mixes many interactions. Re-run profvis while doing exactly one thing (e.g., select one SOM node) to separate per-output costs from framework noise.

Summary

The app is spending most of its time converting ggplots to plotly and redrawing plots that do not need to update. The biggest wins will come from caching plotly objects, reducing the number of plotly outputs, and stopping hidden/irrelevant outputs from firing. After that, stabilizing dimension-reduction caches and precomputing joined/reshaped data will further reduce per-render cost.
