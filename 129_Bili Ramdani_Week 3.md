# *Case Analysis: Modul 1 : Kasus Penyakit Dengue dan Konsep Desain Vaksin (GSE18090)*
### 1. PENDAHULUAN
#### Latar Belakang
Vaksin influenza trivalen merupakan strategi utama dalam pencegahan infeksi influenza dan diketahui tidak hanya meningkatkan pembentukan antibodi, tetapi juga memicu perubahan ekspresi gen pada sel darah perifer sebagai bagian dari respons imun. Analisis transkriptomik memungkinkan evaluasi dinamika perubahan ekspresi gen secara global pada beberapa titik waktu setelah vaksinasi. Dataset GSE48018 menyediakan data ekspresi gen dari pria sehat sebelum serta pada hari ke-1, ke-3, dan ke-14 setelah vaksinasi influenza, sehingga dapat digunakan untuk mengidentifikasi gen yang mengalami upregulation dan downregulation, menentukan 50 Differentially Expressed Genes (DEGs) teratas, serta menganalisis enrichment Gene Ontology (GO) dan KEGG Pathway guna memahami mekanisme molekuler respons vaksin.
#### Tujuan
Membandingkan tingkat ekspresi gen antara dua kondisi biologis (sebelum vs sesudah vaksinasi)  

### 2. Metode
- Dataset: GSE48018 (Vaccine vs Baseline)
- Platform: Microarray Illumina HumanHT-12 V4.0.
- Analisis Statistik: Menggunakan pendekatan statistik limma (*Linear Models for Microarray Data*)
- Preprocessing: Kontrol kualitas data, transformasi log2, dan inspeksi normalisasi.
- Visualisasi Data:
  - Volcano Plot: Melihat signifikansi statistik vs fold change.
  - Heatmap: Pola klastering 50 gen teratas.
  - UMAP: Reduksi dimensi untuk melihat sebaran sampel.
- Tahapan Analisis dengan Perangkat Lunak "R"  
  Tahap 1 - Persiapan Lingkungan Kerja  
  #1. Menginstall BiocManager

      if (!require("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
      }
  #2. Menginstall paket Bioconductor (GEOquery & limma)
  Mengambi data dari database GEO/GEOquery

      BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = 
      FALSE)

      BiocManager::install("illuminaHumanv4.db", ask = FALSE, update = 
      FALSE)
  #3. Menginstall paket CRAN untuk visualisasi dan manipulasi data
      phetmap: heatmap ekspresi gen  
      ggplot2: grafik (volcano plot)  
      dplyr: manipulasi tabel data
      umap: grafik (plot UMAP)
  
      install.packages(c("pheatmap", "ggplot2", "dplyr"))
      if (!requireNamespace("umap", quietly = TRUE)) { 
      install.packages("umap") 
      }
  #4. Memanggil library

      library(GEOquery) 
      library(limma) 
      library(pheatmap) 
      library(ggplot2) 
      library(dplyr) 
      library(illuminaHumanv4.db) 
      library(AnnotationDbi) 
      library(umap)

  Tahap 2 - PENGAMBILAN DATA DARI GEO

      gset <- getGEO("GSE48018", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

  Tahap 3 - PRE-PROCESSING DATA EKSPRESI

      ex <- exprs(gset)
  - Melakukan log2 transformasi

        qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = 
        TRUE))
        LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
        if (LogTransform) {
          ex[ex <= 0] <- NA #Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA 
          ex <- log2(ex) 
          }
  Tahap 4 - MENDEFINISIKAN KELOMPOK SAMPEL
  Mengambil kolom ‘treatment:ch1’ yang berisi info vaksin vs baseline  

      group_info <- pData(gset)[["treatment:ch1"]]
      groups <- make.names(group_info) 
  Mengubah data kategorik menjadi faktor

      gset$group <- factor(groups)
  Melihat kategori unik dalam faktor

      nama_grup <- levels(gset$group) 
      print(nama_grup)
  Tahap 5 - DESIGN MATRIX (KERANGKA STATISTIK)
  Membuat matriks desain untuk model linear
  ~0 berarti TANPA intercept (best practice limma)

      design <- model.matrix(~0 + gset$group)
  Memberi nama kolom agar mudah dibaca

      colnames(design) <- levels(gset$group)
  Menentukan perbandingan biologis (vaksinasi vs baseline)


      grup_vaksin <- "trivalent.influenza.vaccination" 
      grup_baseline <- "baseline"


  Tahap 5 - ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
  Membangun model linear untuk setiap gen

      fit <- lmFit(ex, design)
  Mendefinisikan perbandingan antar grup

      Contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels 
      = design)

  Menerapkan kontras ke model

      fit2 <- contrasts.fit(fit,contrast_matrix)
  Empirical Bayes untuk menstabilkan estimasi varians

      fit2 <- eBayes(fit2)
  
  Mengambil hasil akhir DEG

      topTableResults <- topTable( 
      fit2, 
      adjust = "fdr", 
      sort.by = "B", 
      number = Inf, 
      p.value = 0.05 
      ) 
      head(topTableResults)
  
  Tahap 6 - ANOTASI NAMA GEN
  Mengambil ID probe dari hasil DEG

      probe_ids <- rownames(topTableResults)
  Gunakan database illuminaHumanv4

      gene_annotation <- AnnotationDbi::select( 
      illuminaHumanv4.db, 
      keys = probe_ids, 
      columns = c("SYMBOL", "GENENAME"), 
      keytype = "PROBEID"
  Gabungkan dengan hasil limma

      topTableResults$PROBEID <- rownames(topTableResults)
      topTableResults <- merge( 
      topTableResults, 
      gene_annotation, 
      by = "PROBEID", 
      all.x = TRUE 
      )
      head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

   Tahap 7 - VISUALISASI VOLCANO PLOT

      volcano_data <- data.frame( 
      logFC = topTableResults$logFC, 
      adj.P.Val = topTableResults$adj.P.Val, 
      Gene = topTableResults$SYMBOL 
      )
  Klasifikasi status gen

      volcano_data$status <- "NO" 
      volcano_data$status[volcano_data$logFC > 0.1 & volcano_data$adj.P.Val 
      < 0.05] <- "UP" 
      volcano_data$status[volcano_data$logFC < -0.1 & volcano_data$adj.P.Val 
      < 0.05] <- "DOWN"
  Visualisasi

      ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = 
      status)) + 
        geom_point(alpha = 0.6) + 
        scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = 
      "red")) + 
        geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") + 
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
        theme_minimal() + 
        ggtitle("Volcano Plot DEG Vaksin Flu") 
                      
  Tahap 8 - VISUALISASI HEATMAP
  50 gen paling signifikan berdasarkan adj.P.Val

 
      topTableResults <- topTableResults[ 
        order(topTableResults$adj.P.Val), 
      ] 
      top50 <- head(topTableResults, 50)
      mat_heatmap <- ex[top50$PROBEID, ]
  
      gene_label <- ifelse( 
        is.na(top50$SYMBOL) | top50$SYMBOL == "",
        top50$PROBEID,
        top50$SYMBOL
      )

      rownames(mat_heatmap) <- gene_label

  Pembersihan data (WAJIB agar tidak error hclust) 

      mat_heatmap <- mat_heatmap[ 
        rowSums(is.na(mat_heatmap)) == 0, 
      ]
  Menghapus gen dengan varians nol

      gene_variance <- apply(mat_heatmap, 1, var) 
      mat_heatmap <- mat_heatmap[gene_variance > 0, ]
      annotation_col <- data.frame( 
        Group = gset$group 
      ) 
      rownames(annotation_col) <- colnames(mat_heatmap)
  Visualisasi heatmap

      pheatmap( 
        mat_heatmap, 
        scale = "row",        
        annotation_col = annotation_col, 
        show_colnames = FALSE,         
        show_rownames = TRUE, 
        fontsize_row = 7,
        clustering_distance_rows = "euclidean", 
        clustering_distance_cols = "euclidean", 
        clustering_method = "complete", 
        main = "Top 50 Differentially Expressed Genes" 
        )
  Tahap 9 - VISUALISASI UMAP

      umap_input <- t(ex)
      umap_df <- data.frame( 
        UMAP1 = umap_result$layout[, 1], 
        UMAP2 = umap_result$layout[, 2], 
        Group = gset$group 
      )

      ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) + 
        geom_point(size = 3, alpha = 0.8) + 
        theme_minimal() + 
        labs( 
          title = "UMAP Plot: Baseline vs Vaccinated", 
          x = "UMAP 1", 
          y = "UMAP 2" 
        ) 

   Tahap 10 - MENYIMPAN HASIL
    Menyimpan hasil analisis ke file CSV
  
      write.csv(topTableResults, "Hasil_GSE48018_DEG.csv") 
      message("Analisis selesai. File hasil telah disimpan.")
            
### 3. HASIL
Hasil analisis tertera pada file berikut: 
  
      
