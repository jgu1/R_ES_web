$(function() {
  createGraph();
});

function createGraph() {
  
  //var x_axis_scale = d3.scale.ordinal().rangeBands([0, w])
 
  //var xscale = d3.scale.linear().range([0, w]);
  //var yscale = d3.scale.linear().range([h, 0]);

  var sub_clusterscallback = function(data){
    var gene_p_qs = data.gene_p_qs
    var serisables = data.serisables
    var gene_descriptions = data.gene_descriptions;  
 
    for (var i = 0; i < serisables.length; i++){
        var serisable = serisables[i];
        var row_comb = serisable[0];
        var cols = serisable[1];
        var num_seeds = serisable[2];       
 
        var sub_cluster_gene_p_qs = {};
        var sub_cluster_sorted_pair_names = new Array();
        var sub_cluster_filtered_gene_names = new Array(cols.length);
        var sub_cluster_gene_descriptions = new Array(cols.length);

        // loop through all rows
        for (var i_pairname = 0; i_pairname < row_comb.length; i_pairname++){
            var curr_pairname = row_comb[i_pairname];
            sub_cluster_sorted_pair_names.push(curr_pairname)
            var curr_all_genes = gene_p_qs[curr_pairname];
            curr_pair_genes_arr = new Array();
            // loop though all columns
            for (var i_pass_idx = 0; i_pass_idx < cols.length; i_pass_idx++){
                var curr_pass_idx = cols[i_pass_idx];
                var curr_pass_gene = curr_all_genes[curr_pass_idx];
                curr_pair_genes_arr.push(curr_pass_gene);
                // populate gene names in the first row
                if (curr_pass_gene[2] != "dummy_gene"){
                    sub_cluster_filtered_gene_names[i_pass_idx] = curr_pass_gene[2];
                }
                sub_cluster_gene_descriptions[i_pass_idx] = gene_descriptions[curr_pass_idx]; 
            }
            sub_cluster_gene_p_qs[curr_pairname] = curr_pair_genes_arr;
        }
          
        draw_subcluster(sub_cluster_gene_p_qs,sub_cluster_filtered_gene_names, sub_cluster_sorted_pair_names,sub_cluster_gene_descriptions,num_seeds); 
    }
    // hide the discover_sub_clusters button  
    //d3.select("#sub_clusters_button").style("display","none");
  }
  
  function draw_subcluster(gene_p_qs,geneNames,pairNames,geneDescriptions,num_seeds){
    var matrix = [];
    var size = 15;
    var selectedTextColor = "green"
    pairNames.sort();
    for (var i = 0; i<pairNames.length;i++) {
        pairName = pairNames[i];
        var genes= gene_p_qs[pairName]
        var n = genes.length
        curr_pair = []
        genes.forEach(function(gene){
            curr_pair.push(gene);
        });          
        matrix.push(curr_pair)
    
    }
    var max_pval = 5;
    var min_pval = 0;
    var color_scale=d3.scaleLinear()
                        .domain([min_pval,max_pval])
                        .range(['#ffeda0','#f03b20']);

    var height_pair = matrix.length * size;
    var width_pair= matrix[0].length * size;
    var margin = 120;
 
    var sub_clusters = d3.select("#sub_clusters");
    //sub_clusters.html(""); 
    //var chartHeight = chart.style("height");
    //var chartHeight = parseInt(chartHeight.substring(0,chartHeight.length-2));
    if (height_pair +  margin > 500){
        sub_clusters.attr("height",500);
    }else{
        sub_clusters.attr("height", height_pair + margin);
    }
  
    //var w = width_pair - 2 * margin, h = height_pair - 2 * margin;
    var w = width_pair, h = height_pair;
    //var svg = d3.select("#chart")
    var wrapper_div = sub_clusters.append("div")
                      .attr("width",1000)
                      .style("overflow","scroll");

    var svg = wrapper_div.append("svg")
                    .attr("width", w + 5 * margin)
                    .attr("height",height_pair + margin)
                    //.style("overflow","scroll")
                  .append("svg:g")
                        .attr("transform", "translate(" + 5 * margin + ", " + margin + ")");

    svg.append("rect")
      .attr("class", "background")
      .attr("width", width_pair)
      .attr("height", height_pair);

    var pair = svg.selectAll(".pair")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair")
      .attr("transform", function(d, i) { return "translate(0," + i*size + ")"; })
      .each(pair);

    pair.append("line")
        .attr("x2",w)
        .style("stroke","#fff");
 
    pair.append("text")
        .attr("x",-6)
        .attr("y",7)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .text(function(d,i) { return pairNames[i]; })
        .on("click", function(d,i){ 
            sortAlongRow(i);
        });

    svg.selectAll(".column").data([]).exit().remove();
    var column =svg.selectAll(".column")
              .data(geneNames)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*size + ")rotate(-90)"; });
               
    column.append("line")
        .attr("x1",-height_pair)
        .style("stroke","#fff");

    var columnText = column.append("text")
      .attr("x", 6)
      .attr("y", 7)
      .attr("width",size)
      .attr("height",size)
      .attr("dy", ".32em")
      .text(function(d, i) { return geneNames[i]; })
      .on("click", function(d,i) {sortAlongColumn(i);}) 
      .on("contextmenu",function(d,i){
            d3.event.preventDefault();
            var currColor = this.getAttribute("fill");
            if (currColor != selectedTextColor){
                this.setAttribute("fill",selectedTextColor);
            }else{
                this.setAttribute("fill",null);
            }
        });

    columnText
      .append("title")
        .text(function(d,i){return geneNames[i] +": " + geneDescriptions[i]; });
 
  function pair(pair) {
    var cell = d3.select(this).selectAll(".gene")
        .data(pair)
      .enter().append("rect")
        .attr("class", "gene")
        .attr("x", function(d,i) { return i*size; })
        .attr("width", size)
        .attr("height", size)
        .attr("title",function(d){return d[3]})
        .style("fill", 
            function(d) { 
                if(parseFloat(d[3]) < 0){
                    return '#000'
                }else{
                    return color_scale(-Math.log10(parseFloat(d[3]))) 
                }
            })
     .on("click",function(d){
            //d3.json("/detail?GWAS="+d[0]+"&eQTL="+d[1]+"&gene="+d[2],detailcallback);
            d3.select("#detail").remove();
            wrapper_div.append("div")
               .attr("id","detail")
            d3.json("/detail?gene="+d[2]+"&pairNames="+pairNames,detailcallback);
            //plotDetail(d);
            }
           ) 
      .append("title")
            .text(function(d){
                return "pval:" + d[3] + "\nqval:" + d[4] + "\nGWAS:" + d[0] + "\neQTL:" + d[1] + "\ngene:" + d[2];
            });
       
   
  }

    function getColumnINegLog10(i){
      var index=0;
      var negLog10=[];
      var cd;
      for(index=0; index<matrix.length; index++){
           SNPs=matrix[index];
           if(SNPs.length>i){
                if (parseFloat(SNPs[i][3]) < 0){
                    negLog10.push(-1);
                }else{
                     negLog10.push(-Math.log10(parseFloat(SNPs[i][3])));
                }
           }else{
                negLog10.push(0);
           }
      }
      return negLog10;
    }

    function sortAlongColumn(i){
       var t = svg.transition().duration(3000);
       var ColumnIData=getColumnINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(matrix.length).sort(function(a,b){   //index array pointing to real data arra
            return ColumnIData[b] - ColumnIData[a];
       });
       t.selectAll(".pair")
        .attr("transform", function(d,i) { return "translate(0,"+sorted.indexOf(i) * size+")";});
    }

 
    function getRowINegLog10(i){
      var index=0;
      var negLog10=[];
      var rowIRaw = matrix[i];
      for (var i=0; i<rowIRaw.length; i++){
        gene = rowIRaw[i];
        if (parseFloat(gene[3]) < 0){
            negLog10.push(-1);    
        }else{
            negLog10.push(-Math.log10(parseFloat(gene[3])));  
        }
      }
      return negLog10;
    }


    function sortAlongRow(i){
       var t = svg.transition().duration(3000);
       var rowIData=getRowINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(rowIData.length).sort(function(a,b){   //index array pointing to real data array
            return rowIData[b] - rowIData[a];
       });
       // sort column names
       t.selectAll(".column")
        .attr("transform", function(d,i) { return "translate(" + sorted.indexOf(i) * size+")rotate(-90)";});
       // sort rects in each row 
       var SNPs = t.selectAll(".pair")
                    .selectAll(".gene")
                    .attr("x", function(d,i) { return sorted.indexOf(i)*size; });
    }  
    var Manhattan_btn = wrapper_div.append('button')
                        .attr("type","button")
                        .text("draw Manhattan Plots")
                        .style("margin-left","600px")
                        .on("click",function(d){

                            d3.select("#Manhattan").remove();
                            var Manhattan = wrapper_div.append("div")
                                .attr("id","Manhattan")


                            var selected_geneNames = columnText._groups[0].map(function(text){
                                var oH = text.outerHTML;
                                var selectedTag = 'fill="' + selectedTextColor+'"'; 
                                if (oH.indexOf(selectedTag)>-1){
                                    var geneName = text.__data__;
                                    return geneName;
                                }
                            });
                            Manhattan.attr("geneNames",selected_geneNames)
                            Manhattan.attr("pairNames",pairNames)

                            d3.json("/Manhattan?geneNames="+selected_geneNames+"&pairNames="+pairNames,Manhattancallback);
                            }
                        ) ;
    var a = 1;  
  }

  function append_Manhattan_customizer(parentNode){
    d3.select("#show_unaligned_SNPs").remove();
    d3.select("#txt_show_unaligned_SNPs").remove();
    d3.select("#show_untagged_SNPs") .remove();
    d3.select("#txt_show_untagged_SNPs").remove();
    d3.select("#show_all_genes").remove();
    d3.select("#txt_show_all_genes").remove();
    d3.select("#txt_GSNP_cutoff").remove();
    d3.select("#GSNP_cutoff").remove();
 
    var alignedValue = function(d){return d[4]};
    var taggedValue  = function(d){return d[5]};
    var dotGWAS_size = 3.5;
    var recteQTL_size = 5;
    var invisible_size = 0.1;

    parentNode.append("text")
         .text("GSNP_cutoff")
         .attr("id","txt_GSNP_cutoff")
         .style("margin-left","600px");
    parentNode.append("input")
         .attr("type","text")
         .attr("size","10")
         .attr("id","GSNP_cutoff");


    parentNode.append("text")
         .text("unaligned")
         .attr("id","txt_show_unaligned_SNPs");
    var checkbox_unaligned= parentNode.append("input")
         .attr("type","checkbox")
         .attr("checked",null)
         .attr("id","show_unaligned_SNPs")
         .on("change",function(d){
            d3.selectAll(".dotGWAS") .data([]).exit().remove();
            d3.selectAll(".recteQTL").data([]).exit().remove();

            var zoom_domain_min = parseInt(parentNode.attr("zoom_domain_min"));
            var zoom_domain_max = parseInt(parentNode.attr("zoom_domain_max"));
            var zoom_range_min  = parseInt(parentNode.attr("zoom_range_min"));
            var zoom_range_max  = parseInt(parentNode.attr("zoom_range_max"));
            var Manhattan_width = parseInt(parentNode.attr("Manhattan_width"));
            var Manhattan_height= parseInt(parentNode.attr("Manhattan_height"));
            var geneNames       = parentNode.attr("geneNames");
            var pairNames       = parentNode.attr("pairNames");
            var GSNP_cutoff = d3.select("#GSNP_cutoff").property("value");                
            if ( ! GSNP_cutoff.trim()){
                GSNP_cutoff = 0.001;    
            }

            d3.json("/Manhattan_appendSNPs?zoom_domain_min="+zoom_domain_min+"&zoom_domain_max="+zoom_domain_max+"&zoom_range_min="+zoom_range_min+"&zoom_range_max="+zoom_range_max+"&Manhattan_height="+Manhattan_height+"&geneNames="+geneNames+"&pairNames="+pairNames+"&GSNP_cutoff="+GSNP_cutoff,Manhattan_appendSNPs_callback);
            
         });

    parentNode.append("text")
         .text("untagged")
         .attr("id","txt_show_untagged_SNPs");
    var checkbox_untagged = parentNode.append("input")
         .attr("type","checkbox")
         .attr("checked",null)
         .attr("id","show_untagged_SNPs")
         .on("change",function(d){
            d3.selectAll(".dotGWAS") .data([]).exit().remove();
            d3.selectAll(".recteQTL").data([]).exit().remove();

            var zoom_domain_min = parseInt(parentNode.attr("zoom_domain_min"));
            var zoom_domain_max = parseInt(parentNode.attr("zoom_domain_max"));
            var zoom_range_min  = parseInt(parentNode.attr("zoom_range_min"));
            var zoom_range_max  = parseInt(parentNode.attr("zoom_range_max"));
            var Manhattan_width = parseInt(parentNode.attr("Manhattan_width"));
            var Manhattan_height= parseInt(parentNode.attr("Manhattan_height"));
            var geneNames       = parentNode.attr("geneNames");
            var pairNames       = parentNode.attr("pairNames");
            var GSNP_cutoff = d3.select("#GSNP_cutoff").property("value");                
            if ( ! GSNP_cutoff.trim()){
                GSNP_cutoff = 0.001;    
            }
            d3.json("/Manhattan_appendSNPs?zoom_domain_min="+zoom_domain_min+"&zoom_domain_max="+zoom_domain_max+"&zoom_range_min="+zoom_range_min+"&zoom_range_max="+zoom_range_max+"&Manhattan_height="+Manhattan_height+"&geneNames="+geneNames+"&pairNames="+pairNames+"&GSNP_cutoff="+GSNP_cutoff,Manhattan_appendSNPs_callback);
            
         });


    parentNode.append("div")
        .html("<br/>");


    parentNode.append("text")
         .text("show_all_genes")
         .attr("id","txt_show_all_genes")
         .style("margin-left","600px");
    var checkbox_show_all_genes = parentNode.append("input")
        .attr("type","checkbox")
        .attr("checked",null)
        .attr("id","show_all_genes")
        .on("change",function(d){
        if (this.checked){
            var zoom_domain_min = parseInt(parentNode.attr("zoom_domain_min"));
            var zoom_domain_max = parseInt(parentNode.attr("zoom_domain_max"));
            var zoom_range_min  = parseInt(parentNode.attr("zoom_range_min"));
            var zoom_range_max  = parseInt(parentNode.attr("zoom_range_max"));
            var Manhattan_height= parseInt(parentNode.attr("Manhattan_height"));


            d3.json("/Manhattan_appendGenes?zoom_domain_min="+zoom_domain_min+"&zoom_domain_max="+zoom_domain_max+"&zoom_range_min="+zoom_range_min+"&zoom_range_max="+zoom_range_max+"&Manhattan_height="+Manhattan_height,Manhattan_appendGenes_callback);

     
            //var closest_gene =svg.append("line")
                    
        }else{
            d3.selectAll(".gene_within_domain").remove(); 
        }
         });
  
    checkbox_unaligned.property('checked', false);
    checkbox_untagged.property('checked', false);
    checkbox_show_all_genes.property('checked', false);
 
    parentNode.append("div")
        .html("<br/>");
  }


  var detailcallback = function(data){
    var matrix = [];
    var pairNames = [];
    var size = 15; 
    var gene = data.gene;
    var all_SNPs_list = data.all_SNPs_list
    var pair_SNP_dict = data.pair_SNP_dict;
    var curr_pair = [];

    for (var pairName in pair_SNP_dict) {
        pairNames.push(pairName);
        var SNPs= pair_SNP_dict[pairName];
        var n = SNPs.length;
        curr_pair = []
        SNPs.forEach(function(SNP){
            curr_pair.push(SNP);
        });  
        matrix.push(curr_pair); 
    }
    var detail_height = matrix.length * size;
    var detail_width = matrix[0].length * size; // should not be null, because at least one gene has SNPs
    var margin = {top: 120, right: 200, bottom: 10, left: 600};
    
    d3.select("#detail").html("");
    d3.select("#detail").append("div") 
        .attr("id","gene_of_interest")
        .text("Gene of Interest: " + gene)
        .style("margin-left",margin.left +"px");
   
    var detail = d3.select("#detail");

    var canvas2=detail
           .append("svg")
           .attr("width", detail_width + margin.left)
           .attr("height", detail_height + margin.top)
   
    //d3.select("#gene_of_interest")
    
    //canvas2.html("");// clear canvas2 before drawing a new one
    var svg2=canvas2
       .append("g")
       .attr("transform", "translate(" + margin.left + "," + margin.top + ")"); 

    var max_pval = 5;
    var min_pval = 0;
    var color_scale=d3.scaleLinear()
                        .domain([min_pval,max_pval])
                        .range(['#ffeda0','#f03b20']);

    var height_detail = matrix.length * size;
    var width_detail= matrix[0].length * size;

    if (height_detail +  margin.top > 500){
        detail.attr("height",500);
    }else{
        detail.attr("height", height_detail + margin.top);
    }

    svg2.append("rect")
      .attr("class", "background")
      .attr("width", width_detail)
      .attr("height", height_detail);

    var pair = svg2.selectAll(".pair2")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair2")
      .attr("transform", function(d, i) { return "translate(0," + i*size + ")"; })
      .each(pair2);

     pair.append("line")
        .attr("x2",width_detail)
        .style("stroke","#fff");

     pair.append("text")
        .attr("x",-6)
        .attr("y",7)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .text(function(d,i) { return pairNames[i]; })
        .on("click", function(d,i){ 
            sortAlongRow(i);
        });

    svg2.selectAll(".column").data([]).exit().remove();
    var column =svg2.selectAll(".column")
              .data(all_SNPs_list)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*size + ")rotate(-90)"; });
    
      
    column.append("line")
        .attr("x1",-height_detail)
        .style("stroke","#fff");
    

    column.append("text")
      .attr("x", 6)
      .attr("y", 7)
      .attr("width",size)
      .attr("height",size)
      .attr("dy", ".32em")
      .text(function(d, i) { return all_SNPs_list[i]; }) 
      .on("click", function(d,i) {sortAlongColumn(i);}) ; 
   


 
    function pair2(pair) {
        var cell = d3.select(this).selectAll(".SNP")
            .data(pair)
          .enter().append("rect")
            .attr("class", "SNP")
            .attr("x", function(d,i) { return i*size; })
            .attr("width", size)
            .attr("height", size)
            .style("fill", 
                function(d) { 
                    if(parseFloat(d[2]) < 0){
                        return '#000'
                    }else{
                        return color_scale(-Math.log10(parseFloat(d[2]))) 
                    }
                })
            .append("title").text(function(d){
             if (parseFloat(d[2]) < 0){
                return '';   
             }else{
                return "GWAS pval: " + d[2] + "\neQTL pval: " + d[3] + "\nGWAS-SNP name: " + d[0] + "\neQTL-SNP name: " + d[1] ;
                //return "GWAS pval :" + d[2] + "\nGWAS_SNP:" + d[0] + "\neQTL_SNP :" + d[1] ;
             }
            });  
    }
 
    function getRowINegLog10(i){
      var index=0;
      var negLog10=[];
      var rowIRaw = matrix[i];
      for (var i=0; i<rowIRaw.length; i++){
        SNP = rowIRaw[i];
        if (parseFloat(SNP[2]) < 0){
            negLog10.push(-1);    
        }else{
            negLog10.push(-Math.log10(parseFloat(SNP[2])));  
        }
      }
      return negLog10;
    }


    function sortAlongRow(i){
       var t = svg2.transition().duration(3000);
       var rowIData=getRowINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(rowIData.length).sort(function(a,b){   //index array pointing to real data array
            return rowIData[b] - rowIData[a];
       });
       // sort column names
       t.selectAll(".column")
        .attr("transform", function(d,i) { return "translate(" + sorted.indexOf(i) * size+")rotate(-90)";});
       // sort rects in each row 
       var SNPs = t.selectAll(".pair2")
                    .selectAll(".SNP")
                    .attr("x", function(d,i) { return sorted.indexOf(i)*size; });
    }  

    function getColumnINegLog10(i){
      var index=0;
      var negLog10=[];
      var cd;
      for(index=0; index<matrix.length; index++){
           SNPs=matrix[index];
           if(SNPs.length>i){
                if (parseFloat(SNPs[i][2]) < 0){
                    negLog10.push(-1);
                }else{
                     negLog10.push(-Math.log10(parseFloat(SNPs[i][2])));
                }
           }else{
                negLog10.push(0);
           }
      }
      return negLog10;
    }

    function sortAlongColumn(i){
       var t = svg2.transition().duration(3000);
       var ColumnIData=getColumnINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(matrix.length).sort(function(a,b){   //index array pointing to real data arra
            return ColumnIData[b] - ColumnIData[a];
       });
       t.selectAll(".pair2")
        .attr("transform", function(d,i) { return "translate(0,"+sorted.indexOf(i) * size+")";});
    }  
  }


  var Manhattancallback = function(data){
    //var location_pval_chrom_SNPlist_dict = data.location_pval_chrom_SNPlist_dict; 
    var GWAS_SNPlist_dict                = data.GWAS_SNPlist_dict
    var eQTL_gene_SNPlist_dict           = data.eQTL_gene_SNPlist_dict
    var Manhattan_pairNames              = data.Manhattan_pairNames;    
    var Manhattan_geneNames              = data.Manhattan_geneNames;
    var chrom_starts_data                = data.chrom_starts;
    var gene_location_dict               = data.gene_location_dict;  
    //var eQTL_SNPlist_dict                = data.eQTL_SNPlist_dict;

    var margin = {top: 50, right: 20, bottom: 30, left: 40, between: 50};
    var total_width  = 1000;
    var total_height = 180; 
    var Manhattan_width = total_width - margin.left - margin.right; // should not be null, because at least one gene has SNPs
    var Manhattan_height = total_height - margin.top  - margin.bottom;
    var xScaleMax = 3500000000;
    var yScaleMax = 8;
    var yScaleMax_eQTL = 16;
    var dotGWAS_size = 3.5;
    var recteQTL_size = 5;


 
    d3.select("#Manhattan").html("");
    var Manhattan = d3.select("#Manhattan");
    Manhattan.attr("Manhattan_width",Manhattan_width);
    Manhattan.attr("Manhattan_height",Manhattan_height);
    Manhattan.attr("zoom_domain_min",0);
    Manhattan.attr("zoom_domain_max",xScaleMax);
    Manhattan.attr("zoom_range_min",0);
    Manhattan.attr("zoom_range_max",Manhattan_width);

    append_Manhattan_customizer(Manhattan);
   
    // setup x 
    var xValue = function(d) { return d[2];}, // data -> value
        xScale = d3.scaleLinear().range([0, Manhattan_width]), // value -> display
        xMap = function(d) { return xScale(xValue(d));}, // data -> display
        xAxis = d3.axisBottom().scale(xScale);

    // setup y
    var yValue = function(d) {  var pval = parseFloat(d[3]); 
                                return -Math.log10(pval);
                             }, // data -> value
        yScale = d3.scaleLinear().range([Manhattan_height, 0]), // value -> display
        yMap = function(d) {
            var yValue_direct = yValue(d);
            var yScale_input  = yValue_direct;
            if (yScale_input > yScaleMax){  // if pval is too large, it will go beyond plot, truncate it
                yScale_input = yScaleMax;
            }
            return yScale(yScale_input);}, // data -> display
        yAxis = d3.axisLeft().scale(yScale);

    // setup x for eQTL

    xMap_eQTL = function(d) {
      var xMap_part =  xScale(xValue(d));
      var shift     = recteQTL_size;
      return xMap_part - shift;
    }; 

    // setup y for eQTL
    var yScale_eQTL = d3.scaleLinear().range([Manhattan_height, 0]), // value -> display
        yMap_eQTL = function(d) {
            var yValue_direct = yValue(d);
            var yScale_input  = yValue_direct;
            if (yScale_input > yScaleMax_eQTL){  // if pval is too large, it will go beyond plot, truncate it
                yScale_input = yScaleMax_eQTL;
            }
            return yScale_eQTL(yScale_input);}; // data -> display
    var yAxis_eQTL = d3.axisLeft().scale(yScale_eQTL);
    //var yAxis_eQTL = d3.axisRight().scale(yScale_eQTL);



    // setup fill color
    var cValue = function(d) {return d[6]}; // color by gene    
    var color = d3.scaleOrdinal(d3.schemeCategory20);

    var alignedValue = function(d){return d[4]};
    var taggedValue  = function(d){return d[5]};

    var tooltip = Manhattan.append("div")
        .attr("class", "tooltip")
        .style("height",100)
        .style("width",200)
        .style("opacity", 0);

    xScale.domain([0,xScaleMax]);
    yScale.domain([0,yScaleMax]);
    yScale_eQTL.domain([0,yScaleMax_eQTL]);

    var xMin = -1,xMax = -1,yMin = -1,yMax = -1;
    Manhattan_pairNames.sort();
    NGENES = Manhattan_geneNames.length;
   
        var gx1 = d3.scaleLinear()
            .range([0,Manhattan_width]);

        var gx2 = d3.scaleLinear()
            .range([0,Manhattan_width]);

        gx1.domain([0, xScaleMax]).nice();
        gx2.domain([0, xScaleMax]).nice();

        var zoomX = d3.zoom().scaleExtent([0.2, 2000000]).on("zoom", zoomedX);
        //Zoom in v4


        function zoomedX(){

            var tr = d3.event.transform;
            gx2 = tr.rescaleX(gx1);
            
            var zoom_domain_min = gx2.domain()[0];
            var zoom_domain_max = gx2.domain()[1];
            var zoom_range_min  = gx2.range()[0];
            var zoom_range_max  = gx2.range()[1];
            Manhattan.attr("zoom_domain_min",zoom_domain_min);
            Manhattan.attr("zoom_domain_max",zoom_domain_max);
            Manhattan.attr("zoom_range_min",zoom_range_min);
            Manhattan.attr("zoom_range_max",zoom_range_max); 

            //console.log("gx2: x,y:"+x+ ' '+ y + ' ')
            //console.log("###gx1: x,y:"+gx1.domain()[0]+ ' '+ gx1.domain()[1])
            Manhattan.selectAll(".Manhattan_group").selectAll(".dotGWAS")
                .attr("cx",xcoord);

            Manhattan.selectAll(".Manhattan_group").selectAll(".recteQTL")
                .attr("x",xcoord);
            Manhattan.selectAll(".Manhattan_group").selectAll(".xAxis")
                .call(xAxis.scale(gx2));

            Manhattan.selectAll(".Manhattan_group").selectAll(".chrom_starts")
                .attr("x1",function(d){return gx2(d)})
                .attr("x2",function(d){return gx2(d)});

            Manhattan.selectAll(".Manhattan_group").selectAll(".gene_starts")
                .attr("x",function(d){
                    var location_tuple = gene_location_dict[d];
                    return gx2(xValue(location_tuple))
                })
                .attr("width",function(d){
                  //location_tuple = (gene,chrom,chrom_abs_start,chrom_abs_end,chromStart,chromEnd)
                  var location_tuple = gene_location_dict[d];
                  var width = 1;
                  chrom_abs_start = location_tuple[2]
                  chrom_abs_end   = location_tuple[3]
                  if (chrom_abs_start != null && chrom_abs_end != null){
                    width = gx2(chrom_abs_end) - gx2(chrom_abs_start);
                  }
                  if (width < 1){ width = 1;}
                  return width;
                });
                /*
                .attr("x",function(d){
                  var chromStart_chromEnd = gene_location_dict[d]; 
                  var chromStart = chromStart_chromEnd[0];
                  var x = 0;
                  if (chromStart != null){
                    x = gx2(chromStart);
                  }
                  return x;
                })
               .attr("width",function(d){
                 var chromStart_chromEnd = gene_location_dict[d]; 
                 var chromStart = chromStart_chromEnd[0];
                 var chromEnd   = chromStart_chromEnd[1];
                 var width = 1;
                 if (chromStart != null && chromEnd != null){
                    width = gx2(chromEnd) - gx2(chromStart);
                 }
                 if (width < 1){ width = 1;}
                 return width;
              });
                */
            
            Manhattan.selectAll(".Manhattan_group").selectAll(".gene_within_domain")
                .attr("x1",function(d){
                    var location_tuple = gene_location_dict[d];
                    return gx2(xValue(location_tuple))
                })
                .attr("x2",function(d){
                    var location_tuple = gene_location_dict[d];
                    return gx2(xValue(location_tuple))
                })
                              /*
                .attr("x1",function(d){
                    var chromStart_chromEnd = gene_location_dict[d]; 
                    var chromStart = chromStart_chromEnd[0];
                    return gx2(chromStart);
                    })
                .attr("x2",function(d){
                    var chromStart_chromEnd = gene_location_dict[d]; 
                    var chromStart = chromStart_chromEnd[0];
                    return gx2(chromStart);
                    });
                */

        };
        //Zoom in v4
        var xcoord = function(d){
            return gx2(xValue(d));
        };

        //Zoom in v4


    var Manhattan_svg_width = Manhattan_width + margin.left + margin.right;
    // the height reserved for each GWAS-eQTL pair,including one GWAS and NGENES eQTL
    var Manhattan_svg_height_single_pair = margin.top + margin.bottom + (Manhattan_height + margin.between) * (NGENES + 1);
    var Manhattan_svg_height = Manhattan_svg_height_single_pair * Manhattan_pairNames.length;

    // the centralized svg to catch zooming action
    var Manhattan_svg_container = Manhattan.append("svg")
                                  .attr("width",Manhattan_svg_width)
                                  .attr("height",Manhattan_svg_height)
                                  .call(zoomX)

    var Manhattan_groups = Manhattan_svg_container.selectAll(".Manhattan_group")
                            .data(Manhattan_pairNames)
                            .enter().append("g")
                            .attr("class","Manhattan_group")                          
                            .attr("transform",function(d,i){
                                var shift_for_pair = i * Manhattan_svg_height_single_pair;
                                var shift_tot      = shift_for_pair + margin.top;
                                return "translate(" + margin.left + "," +  shift_tot + ")"
                            });
    Manhattan_groups.append("rect")
                .attr("width", Manhattan_width)
                .attr("height", Manhattan_height)
                .style("fill", "none")
                .style("pointer-events", "all");
    Manhattan_groups.append("g")
          .attr("class","Manhattan_pairName")
         .append("text")
          .text(function(d){return d;});
    
    // a deep copy of Manhattan_geneNames
    var array_GWAS_genes = JSON.parse(JSON.stringify(Manhattan_geneNames));
    array_GWAS_genes.unshift('GWAS');   // add dummy name 'GWAS' at first position 

    Manhattan_groups.selectAll(".chrom_starts").data([]).exit().remove();
    Manhattan_groups.selectAll(".gene_starts").data([]).exit().remove();

    
    var GWAS_or_eQTL_groups = Manhattan_groups.selectAll('.GWAS_or_eQTL_group')
            .data(array_GWAS_genes)
          .enter().append("g")
            .attr("class","GWAS_or_eQTL_group")
            .attr("transform",function(d,i){
                var shift_for_GWAS_or_eQTL = i * (Manhattan_height + margin.between);
                return "translate(0," +  shift_for_GWAS_or_eQTL + ")"
            });

    GWAS_or_eQTL_groups.append("g").selectAll(".chrom_starts")
            .data(chrom_starts_data)
          .enter().append("line")
            .attr("class","chrom_starts")
            .attr("x1",function(d){return xScale(d);})
            .attr("x2",function(d){return xScale(d);})
            .attr("y1",0)
            .attr("y2",Manhattan_height)
            .style("opacity", .2)
            .style("stroke","black");
    
    GWAS_or_eQTL_groups.selectAll(".gene_starts")
            .data(Manhattan_geneNames)
          .enter().append("rect")
            .attr("class","gene_starts")
            .attr("x",function(d){
                var location_tuple = gene_location_dict[d];
                if (location_tuple == null){
                    return 0;
                }
                return xMap(location_tuple)
            })
            .attr("width",function(d){
                //location_tuple = (gene,chrom,chrom_abs_start,chrom_abs_end,chromStart,chromEnd)
                var location_tuple = gene_location_dict[d];
                if (location_tuple == null){
                    return 0;
                }
                var width = 1;
                chrom_abs_start = location_tuple[2]
                chrom_abs_end   = location_tuple[3]
                if (chrom_abs_start != null && chrom_abs_end != null){
                width = gx2(chrom_abs_end) - gx2(chrom_abs_start);
                }
                if (width < 1){ width = 1;}
                return width;
            })
            .attr("y",0)
            .attr("height",Manhattan_height)
            .attr("fill",function(d){
                return color(d);
            });

    GWAS_or_eQTL_groups.append("g")
          .attr("class", "xAxis")
          .attr("transform", "translate(0," + Manhattan_height + ")")
          .call(xAxis)
        .append("text")
          .attr("class", "label")
          .attr("x", Manhattan_width)
          .attr("y", -6)
          .style("text-anchor", "end")
          .text("Chrom loc");
 

    GWAS_or_eQTL_groups.filter(function(d){return d == "GWAS"})
          .append("g")
          .attr("class","yAxis") 
          .call(yAxis);
        
    GWAS_or_eQTL_groups.filter(function(d){return d != "GWAS"})
          .append("g")
          .attr("class","yAxis_eQTL") 
          .call(yAxis_eQTL);

    GWAS_or_eQTL_groups.filter(function(d){return d == "GWAS"})
        .selectAll(".dotGWAS")    
        .data(function(d){
                var parent_Manhattan_group = this.parentNode;
                var pairName = d3.select(this.parentNode).datum();
                return GWAS_SNPlist_dict[pairName];
            })
        .enter().append("circle")
          .attr("class", "dotGWAS")
          .attr("r", dotGWAS_size)
          .attr("cx", xMap)
          .attr("cy", yMap)
          .attr("aligned",alignedValue)
          .attr("tagged",taggedValue)
          .style("fill","black") 
          .on("mouseover", function(d) {
              tooltip.html("");
              tooltip.transition()
                   .duration(200)
                   .style("opacity", .9);
      
              tooltip.style("left", (d3.event.pageX ) + "px")
                   .style("top", function(d){
                        return (d3.event.pageY - 60) + "px"})
                   .html(d[0] + "(" + d[1] + ")"
                        + " <br/>pval: " + yValue(d)
                        + " <br/>aligned: " + alignedValue(d)
                        + " <br/>tagged: " + taggedValue(d)
                        + " <br/>closest gene: " + d[6] 
                        + " <br/>concat_abs: " + d[2]
                        );
              var closest_geneName = d[6];
              var closest_gene_chromStart = d[7];
                
                //var closest_gene =svg.append("line")
                var closest_gene =d3.selectAll(".Manhattan_group").append("line")
                  .attr("class","closest_gene")
                  .attr("x1",gx2(closest_gene_chromStart))
                  .attr("x2",gx2(closest_gene_chromStart))
                  .attr("y1",0)
                  .attr("y2",Manhattan_height) 
                  .style("opacity", 1)
                  .style("stroke","black");
              })
          .on("mouseout", function(d) {
              tooltip.transition()
                   .duration(500)
                   .style("opacity", 0);

              d3.selectAll(".closest_gene").data([]).exit().remove();
              });
     
    GWAS_or_eQTL_groups.filter(function(d){return d != "GWAS"})
        .selectAll(".recteQTL")    
        .data(function(d){
            var parent_Manhattan_group = this.parentNode;
            var pairName = d3.select(this.parentNode).datum();
            var gene_SNPlist_dict = eQTL_gene_SNPlist_dict[pairName];        
            var SNPlist_for_current_gene = gene_SNPlist_dict[d]; 
            return SNPlist_for_current_gene;
        }).enter().append("rect")
          .attr("class", "recteQTL")
          .attr("width", recteQTL_size)
          .attr("height",recteQTL_size)
          .attr("x", xMap_eQTL)
          .attr("y", yMap_eQTL) 
          .attr("aligned",alignedValue)
          .attr("tagged",taggedValue)
          .style("fill", function(d) {
            var gene = cValue(d); 
            return color(gene);                 //gene
            }) 
          .style("fill-opacity",0.9)
          .style("stroke",function(d){return "black"})
          .style("stroke-opacity",0.1)
          .style("stoke-width",function(d){return 5}) 
          .on("mouseover", function(d) {
                  tooltip.html("");
                  tooltip.transition()
                       .duration(200)
                       .style("opacity", .9);
          
                  tooltip.style("left", (d3.event.pageX ) + "px")
                       .style("top", function(d){
                            return (d3.event.pageY - 60) + "px"})
                       .html(d[0] + "(" + d[1] + ")"
                            + " <br/>pval: " + yValue(d) 
                            + " <br/>gene: " + d[6]
                            + " <br/>aligned GSNP: " + d[7]
                            + " <br/>aligned: " + alignedValue(d)
                            + " <br/>tagged: " + taggedValue(d)
                            );
              })
          .on("mouseout", function(d) {
                  tooltip.transition()
                       .duration(500)
                       .style("opacity", 0);
          });

    var legend_block_size = 10;
    var legend = GWAS_or_eQTL_groups.selectAll(".legend")
      .data(color.domain())
    .enter().append("g")
      .attr("class", "legend")
      .attr("transform", function(d, i) { return "translate(0," + i * legend_block_size + ")"; });

    // draw legend colored rectangles
    legend.append("rect")
      .attr("x", Manhattan_width - 10*legend_block_size)
      //.attr("x", 0)
      .attr("width", legend_block_size -1)
      .attr("height", legend_block_size -1)
      .style("fill", color);

    // draw legend text
    legend.append("text")
      .attr("x", Manhattan_width - 9 * legend_block_size)
      .attr("y", legend_block_size -1)
      .style("font-size","9px")
      .style("text-anchor", "start")
      .text(function(d) { return d;})

    toggle_dotGWAS_and_recteQTL();

  }

  var Manhattan_appendGenes_callback = function(data){
    var zoom_domain_min         = data.zoom_domain_min;
    var zoom_domain_max         = data.zoom_domain_max;
    var zoom_range_min          = data.zoom_range_min;    
    var zoom_range_max          = data.zoom_range_max;
    var geneNames               = data.geneNames;
    var Manhattan_height        = data.Manhattan_height;
    var gene_location_dict      = data.gene_location_dict;

    var xValue = function(d) { return d[2];}, // data -> value
        xMap = function(d) { return gx2(xValue(d));} // data -> display


    var gx2 = d3.scaleLinear()
                .range([zoom_range_min,zoom_range_max]);
        gx2.domain([zoom_domain_min, zoom_domain_max]).nice();

    d3.selectAll(".tooptip").data([]).exit().remove();
    var tooltip = d3.select("#Manhattan").append("div")
        .attr("class", "tooltip")
        .style("height",100)
        .style("width",200)
        .style("opacity", 0);


    d3.selectAll(".Manhattan_group").append("g").selectAll(".gene_within_domain")
          .data(geneNames)
          .enter().append("line")
          .attr("class","gene_within_domain")
          .attr("x1",function(d){
                    var location_tuple = gene_location_dict[d];
                    return gx2(xValue(location_tuple))
                })
          .attr("x2",function(d){
                    var location_tuple = gene_location_dict[d];
                    return gx2(xValue(location_tuple))
                })
          /*
          .attr("x1",function(d){
                    var gene_start = gene_location_dict[d];
                    return gx2(gene_location_dict[d]);})
          .attr("x2",function(d){return gx2(gene_location_dict[d]);})
          */
          .attr("y1",0)
          .attr("y2",Manhattan_height / 2) 
          .style("opacity",0.2 )
          .style("stroke","black")
          .on("mouseover", function(d) {
                  tooltip.html("");
                  tooltip.transition()
                       .duration(200)
                       .style("opacity", .9);
          
                  tooltip.style("left", (d3.event.pageX ) + "px")
                       .style("top", function(d){
                            return (d3.event.pageY - 60) + "px"})
                       .html(  " <br/>geneName: " + d
                             + " <br/>chrom: " + gene_location_dict[d][1]
                             + " <br/>chromStart:" + gene_location_dict[d][4]   
                             + " <br/>chromEnd:"   + gene_location_dict[d][5]   
                        );
     
          });


  }

  var Manhattan_appendSNPs_callback = function(data){
    var zoom_domain_min         = data.zoom_domain_min;
    var zoom_domain_max         = data.zoom_domain_max;
    var zoom_range_min          = data.zoom_range_min;    
    var zoom_range_max          = data.zoom_range_max;
    var Manhattan_height        = parseFloat(data.Manhattan_height);
    var GWAS_SNPlist_dict       = data.GWAS_SNPlist_dict;
    var eQTL_gene_SNPlist_dict  = data.eQTL_gene_SNPlist_dict;

    d3.selectAll(".tooptip").data([]).exit().remove();
    var tooltip = d3.select("#Manhattan").append("div")
        .attr("class", "tooltip")
        .style("height",100)
        .style("width",200)
        .style("opacity", 0);

    var color = d3.scaleOrdinal(d3.schemeCategory20);
    var cValue = function(d) {return d[6]}; // color by gene    
    var alignedValue = function(d){return d[4]};
    var taggedValue  = function(d){return d[5]};

    var gx2 = d3.scaleLinear()
                .range([zoom_range_min,zoom_range_max]);
        gx2.domain([zoom_domain_min, zoom_domain_max]).nice();
    var xValue = function(d) { return d[2];}, // data -> value
        xMap = function(d) { return gx2(xValue(d));} // data -> display

    // setup y
    var yScaleMax = 8;
    var dotGWAS_size = 3.5;
    var recteQTL_size = 5;
    var yValue = function(d) {  var pval = parseFloat(d[3]); 
                                return -Math.log10(pval);
                             }, // data -> value
        yScale = d3.scaleLinear().range([Manhattan_height, 0]), // value -> display
        yMap = function(d) {
            var yValue_direct = yValue(d);
            var yScale_input  = yValue_direct;
            if (yScale_input > yScaleMax){  // if pval is too large, it will go beyond plot, truncate it
                yScale_input = yScaleMax;
            }
            return yScale(yScale_input);} // data -> display
    
    yScale.domain([0,yScaleMax]);
    var alignedValue = function(d){return d[4]};
    var taggedValue  = function(d){return d[5]};

    d3.selectAll(".Manhattan_group").selectAll(".GWAS_or_eQTL_group").filter(function(d){return d == "GWAS"})
            .selectAll(".dotGWAS")
          .data(function(d){
                var parent_Manhattan_group = this.parentNode;
                var pairName = d3.select(this.parentNode).datum();
                var SNPlist_for_curr_GWAS = GWAS_SNPlist_dict[pairName];
                return SNPlist_for_curr_GWAS;
            })
          .enter().append("circle")
          .attr("class","dotGWAS")
          .attr("r", dotGWAS_size)
          .attr("cx", xMap)
          .attr("cy", yMap)
          .attr("aligned",alignedValue)
          .attr("tagged",taggedValue)
          .style("fill","black")
          .on("mouseover", function(d) {
              tooltip.html("");
              tooltip.transition()
                   .duration(200)
                   .style("opacity", .9);
      
              tooltip.style("left", (d3.event.pageX ) + "px")
                   .style("top", function(d){
                        return (d3.event.pageY - 60) + "px"})
                   .html(d[0] + "(" + d[1] + ")"
                        + " <br/>pval: " + yValue(d)
                        + " <br/>aligned: " + alignedValue(d)
                        + " <br/>tagged: " + taggedValue(d)
                        + " <br/>closest gene: " + d[6] 
                        );
              var closest_geneName = d[6];
              var closest_gene_chromStart = d[7];
                
                //var closest_gene =svg.append("line")
                var closest_gene =d3.selectAll(".Manhattan_group").append("line")
                  .attr("class","closest_gene")
                  .attr("x1",gx2(closest_gene_chromStart))
                  .attr("x2",gx2(closest_gene_chromStart))
                  .attr("y1",0)
                  .attr("y2",Manhattan_height) 
                  .style("opacity", 1)
                  .style("stroke","black");
              })
          .on("mouseout", function(d) {
              tooltip.transition()
                   .duration(500)
                   .style("opacity", 0);

              d3.selectAll(".closest_gene").data([]).exit().remove();
              });
     

    d3.selectAll(".Manhattan_group").selectAll(".GWAS_or_eQTL_group").filter(function(d){return d != "GWAS"})
            .selectAll(".recteQTL")
          .data(function(d){
                var parent_Manhattan_group = this.parentNode;
                var pairName = d3.select(this.parentNode).datum();
                var gene_SNPlist_dict = eQTL_gene_SNPlist_dict[pairName];
                var SNPlist_for_curr_gene = gene_SNPlist_dict[d];
                return SNPlist_for_curr_gene;
            })
          .enter().append("rect")
          .attr("class", "recteQTL")
          .attr("width", recteQTL_size)
          .attr("height",recteQTL_size)
          .attr("x", xMap)
          .attr("y", yMap) 
          .attr("aligned",alignedValue)
          .attr("tagged",taggedValue)
          .style("fill", function(d) {
            var gene = cValue(d); 
            return color(gene);                 //gene
            }) 
          .style("fill-opacity",0.9)
          .style("stroke",function(d){return "black"})
          .style("stroke-opacity",0.1)
          .style("stoke-width",function(d){return 5}) 
          .on("mouseover", function(d) {
                  tooltip.html("");
                  tooltip.transition()
                       .duration(200)
                       .style("opacity", .9);
          
                  tooltip.style("left", (d3.event.pageX ) + "px")
                       .style("top", function(d){
                            return (d3.event.pageY - 60) + "px"})
                       .html(d[0] + "(" + d[1] + ")"
                            + " <br/>pval: " + yValue(d) 
                            + " <br/>gene: " + d[6]
                            + " <br/>aligned GSNP: " + d[7]
                            + " <br/>aligned: " + alignedValue(d)
                            + " <br/>tagged: " + taggedValue(d)
                            );
              })
          .on("mouseout", function(d) {
                  tooltip.transition()
                       .duration(500)
                       .style("opacity", 0);
          });







    toggle_dotGWAS_and_recteQTL(); 
  }

  function toggle_dotGWAS_and_recteQTL(){
    var dotGWAS_size = 3.5;
    var recteQTL_size = 5;
    var invisible_size = 0.1;

    var alignedValue = function(d){return d[4]};
    var taggedValue  = function(d){return d[5]};


    var show_unaligned = document.getElementById('show_unaligned_SNPs').checked;
    var show_untagged  = document.getElementById('show_untagged_SNPs').checked;

    if (show_unaligned && show_untagged){  
        d3.selectAll(".dotGWAS") .filter(function(d) {return !alignedValue(d);})
            .style("opacity",1)
            .attr("r",dotGWAS_size);
        d3.selectAll(".recteQTL").filter(function(d) {return !alignedValue(d);}) 
            .style("opacity",0.9)
            .attr("width",  recteQTL_size).attr("height", recteQTL_size);
   }else{
        d3.selectAll(".dotGWAS") .filter(function(d) {return !alignedValue(d);})
            .style("opacity",0)
            .attr("r",invisible_size);
        d3.selectAll(".recteQTL").filter(function(d) {return !alignedValue(d);})
            .style("opacity",0)
            .attr("width", invisible_size).attr("height",invisible_size);
   }
   
   if (show_untagged){ 
        d3.selectAll(".recteQTL").filter(function(d) {return !taggedValue(d);})
            .style("opacity",0.9)
            .attr("width", recteQTL_size).attr("height",recteQTL_size);
        d3.selectAll(".dotGWAS") .filter(function(d) {return !taggedValue(d);})
            .style("opacity",1)
            .attr("r",dotGWAS_size);
   
        // there are two types of untagged SNPs: untagged_aligned, untagged_unaligned
        // when showing unaligned_SNPs, both should be brought back to display
        // BUT when "NOT showing unalinged_SNPs", only untagged_aligned should be brought back to display
        if ( ! show_unaligned){
            d3.selectAll(".recteQTL").filter(function(d) {return !alignedValue(d);}) 
                .style("opacity",0)
                .attr("width",invisible_size).attr("height",invisible_size);
            
            d3.selectAll(".dotGWAS").filter(function(d) {return !alignedValue(d);}) 
                .style("opacity",0)
                .attr("r",invisible_size);
        }                
    }else{ 
        d3.selectAll(".recteQTL").filter(function(d) {return !taggedValue(d);})
            .style("opacity",0)
            .attr("width",invisible_size).attr("height",invisible_size);
        d3.selectAll(".dotGWAS") .filter(function(d) {return !taggedValue(d);})
            .style("opacity",0)
            .attr("r",invisible_size);
    }




  }


  function draw_input_fields_ISA(wrapper_div){
    var ISA_input_wrapper_div = wrapper_div.append("div")
                                .attr("id","ISA_input_wrapper_div");

    ISA_input_wrapper_div.append("text")
         .text("filter_ratio");
    ISA_input_wrapper_div.append("input")
         .attr("type","text")
         .attr("size","10")
         .attr("id","ISA_filter_ratio");

    wrapper_div.append("text")
         .text("binarize_cutoff");
    wrapper_div.append("input")
         .attr("type","text")
         .attr("size","10")
         .attr("id","binarize_cutoff");

    wrapper_div.append("text")
         .text("consider_all_genes_in_database")
    wrapper_div.append("input")
         .attr("type","checkbox")
         .attr("checked",null)
         .attr("id","consider_all_genes_in_database");
  }

  function draw_input_fields_PLAID(wrapper_div){
    var PLAID_input_wrapper_div = wrapper_div.append("div")
                                .attr("id","PLAID_input_wrapper_div");
    wrapper_div.append("text")
         .text("binarize_cutoff");
    wrapper_div.append("input")
         .attr("type","text")
         .attr("size","10")
         .attr("id","binarize_cutoff");
   
    wrapper_div.append("text")
         .text("consider_all_genes_in_database")
    wrapper_div.append("input")
         .attr("type","checkbox")
         .attr("checked",null)
         .attr("id","consider_all_genes_in_database");
  }

  var draw_pair = function(draw_pair_json_obj){
    var selectedTextColor = "green"
    var matrix = [];
    var pairNames = [];
    var geneNames = [];
    var size = 15;
    var geneNames = draw_pair_json_obj.filtered_gene_names; 
    var geneDescriptions = draw_pair_json_obj.gene_descriptions;
    var gene_p_qs = draw_pair_json_obj.gene_p_qs; 
    var pairNames = draw_pair_json_obj.sorted_pair_names; 
    for (var i = 0; i<pairNames.length;i++) {
        pairName = pairNames[i];
        var genes= gene_p_qs[pairName]
        var n = genes.length
        curr_pair = []
        genes.forEach(function(gene){
            curr_pair.push(gene);
        });          
        matrix.push(curr_pair)
    
    }   
    var max_pval = 5;
    var min_pval = 0;
    var color_scale=d3.scaleLinear()
                        .domain([min_pval,max_pval])
                        .range(['#ffeda0','#f03b20']);


    var height_pair = matrix.length * size;
    var width_pair= matrix[0].length * size;
    var margin = 120;
    
    var chart = d3.select("#chart");
    //var chartHeight = chart.style("height");
    //var chartHeight = parseInt(chartHeight.substring(0,chartHeight.length-2));
    if (height_pair +  margin > 500){
        chart.attr("height",500);
    }else{
        chart.attr("height", height_pair + margin);
    }
  
    //var w = width_pair - 2 * margin, h = height_pair - 2 * margin;
    var w = width_pair, h = height_pair;
    var chart = d3.select("#chart");
    chart.selectAll("svg").remove(); 
    //var svg = d3.select("#chart")
    var svg = chart.append("svg")
                    .attr("width", w + 5 * margin)
                    .attr("height",height_pair + margin)
                    .style("overflow","scroll")
                  .append("svg:g")
                        .attr("transform", "translate(" + 5 * margin + ", " + margin + ")");



    svg.append("rect")
      .attr("class", "background")
      .attr("width", width_pair)
      .attr("height", height_pair);

   // show_sub_clusters_button = document.getElementById('show_discover_sub_clusters_button').value 
 
    var Manhattan_btn = chart.append('button')
                    .attr("type","button")
                    .text("draw Manhattan Plots")
                    .style("margin-left","600px")
                    .on("click",function(d){
     
                        d3.select("#Manhattan").remove();
                        var Manhattan = chart.append("div")
                            .attr("id","Manhattan")

                        //var GSNP_cutoff = document.getElementById("GSNP_cutoff_global").value;
                    
                        var selected_geneNames = columnText._groups[0].map(function(text){
                            var oH = text.outerHTML;
                            var selectedTag = 'fill="' + selectedTextColor+'"'; 
                            if (oH.indexOf(selectedTag)>-1){
                                var geneName = text.__data__;
                                return geneName;
                            }
                        });
                        Manhattan.attr("geneNames",selected_geneNames)
                        Manhattan.attr("pairNames",pairNames)
                        d3.json("/Manhattan?geneNames="+selected_geneNames+"&pairNames="+pairNames,Manhattancallback);
                        //d3.json("/Manhattan?geneNames="+selected_geneNames+"&pairNames="+pairNames+"&GSNP_cutoff="+GSNP_cutoff,Manhattancallback);
                        }
                    ) ;
              
    var clustering_algs_input_wrapper_div = chart.append("div")
                                            .attr("id","clustering_algs_input_wrapper_div")
                                            .style("margin-left", 5* margin + "px");

    var clustering_algs_wrapper_div = chart.append("div")
                                        .attr("id",'clustering_algs_wrapper_div')
                                        .style("margin-left", 5* margin + "px");
    clustering_algs_wrapper_div.append("text")
         .text("ISA");
    clustering_algs_wrapper_div.append("input")
         .attr("type","checkbox")
         .attr("checked",null)
         .attr("id","clustering_algs_ISA")
         .on("change",function(d){
            clustering_algs_input_wrapper_div.html("");
            if (this.checked) {
                document.getElementById('clustering_algs_PLAID').checked = false;
                draw_input_fields_ISA(clustering_algs_input_wrapper_div);
            }
         });

 
    clustering_algs_wrapper_div.append("text")
         .text("PLAID");
    clustering_algs_wrapper_div.append("input")
         .attr("type","checkbox")
         .attr("checked",null)
         .attr("id","clustering_algs_PLAID")
         .on("change",function(d){
            clustering_algs_input_wrapper_div.html("");
            if (this.checked) {
                document.getElementById('clustering_algs_ISA').checked = false;
                draw_input_fields_PLAID(clustering_algs_input_wrapper_div);
            }

        });
  

                                       
    
//    draw_ISA_input_fields(clustering_algs_input_wrapper_div);

        chart.append("input")
             .attr("type","button")
             .attr("value","Discover sub_clusters")
             .attr("id","sub_clusters_button")
             .style("margin-left",5 * margin +"px")
             .on("click",function(d){

                var sub_clusters = d3.select("#sub_clusters");
                sub_clusters.html("");

                var binarize_cutoff = document.getElementById("binarize_cutoff").value;
                var consider_all_genes_in_database = document.getElementById("consider_all_genes_in_database").checked


                if (document.getElementById("clustering_algs_ISA").checked){
                    var filter_ratio  = document.getElementById("ISA_filter_ratio").value;
                    
                    d3.json("/sub_clusters?alg=ISA&filter_ratio="+filter_ratio + "&binarize_cutoff=" + binarize_cutoff + "&consider_all_genes_in_database="+consider_all_genes_in_database,sub_clusterscallback);
                }else if (document.getElementById("clustering_algs_PLAID").checked){
                    d3.json("/sub_clusters?alg=PLAID&binarize_cutoff=" + binarize_cutoff + "&consider_all_genes_in_database="+consider_all_genes_in_database,sub_clusterscallback);
                }
            });
   // }





    var pair = svg.selectAll(".pair")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair")
      .attr("transform", function(d, i) { return "translate(0," + i*size + ")"; })
      .each(pair);

    pair.append("line")
        .attr("x2",w)
        .style("stroke","#fff");
 
    function getRowINegLog10(i){
      var index=0;
      var negLog10=[];
      var rowIRaw = matrix[i];
      for (var i=0; i<rowIRaw.length; i++){
        gene = rowIRaw[i];
        pval=parseFloat(gene[3]);
        if (pval < 0){
            negLog10.push(-1);
        }else{
            negLog10.push(-Math.log10(pval));  
        }  
      }
      return negLog10;
    }

    /*
    function sortEachRow(row){
        var rects = row.selectAll(".gene")
            .attr("transform",function(d,i){return "translate("+(30-i) * size+")";});
    }
    */

    function sortAlongRow(i){
       var t = svg.transition().duration(3000);
       var rowIData=getRowINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(rowIData.length).sort(function(a,b){   //index array pointing to real data array
        /*
        if(rowSortOrderDesc == true){
            return rowIData[a] - rowIData[b];
        }else{
            return rowIData[b] - rowIData[a];
        }
        */
            return rowIData[b] - rowIData[a];
       });
       // sort column names
       t.selectAll(".column")
        .attr("transform", function(d,i) { return "translate(" + sorted.indexOf(i) * size+")rotate(-90)";});
       // sort rects in each row 
       var genes = t.selectAll(".pair")
                    .selectAll(".gene")
                    .attr("x", function(d,i) { return sorted.indexOf(i)*size; });
    }  

    pair.append("text")
        .attr("x",-6)
        .attr("y",7)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .text(function(d,i) { return pairNames[i]; })
        .on("click", function(d,i){ 
            sortAlongRow(i);
        });
        /*
        .on("click", function(d,i) {
            d3.json("/sortAlongPair?sortAlongPairName="+pairNames[i],draw_pair);
          });
        */
    svg.selectAll(".column").data([]).exit().remove();
    var column =svg.selectAll(".column")
              .data(geneNames)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*size + ")rotate(-90)"; });
               
    column.append("line")
          .style("stroke","#fff");
    
    function getColumnINegLog10(i){
      var index=0;
      var negLog10=[];
      var cd;
      for(index=0; index<matrix.length; index++){
           genes=matrix[index];
           if(genes.length>i){
                columnINumber = parseFloat(genes[i][3]);
                if (columnINumber < 0){
                    negLog10.push(0);
                }else{
                    negLog10.push(-Math.log10(columnINumber));
                }
           }else{
                negLog10.push(0);
           }
      }
      return negLog10;
    }

    function sortAlongColumn(i){
       var t = svg.transition().duration(3000);
       var ColumnIData=getColumnINegLog10(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(matrix.length).sort(function(a,b){   //index array pointing to real data arra
            return ColumnIData[b] - ColumnIData[a];
       });
       t.selectAll(".pair")
        .attr("transform", function(d,i) { return "translate(0,"+sorted.indexOf(i) * size+")";});
  }  


    column.append("line")
        .attr("x1",-w)
        .style("stroke","#fff");


    var columnText = column.append("text")
      .attr("x", 6)
      .attr("y", 7)
      .attr("width",size)
      .attr("height",size)
      .attr("dy", ".32em")
      .text(function(d, i) { return geneNames[i]; }) 
      .on("click", function(d,i) {sortAlongColumn(i);})
      .on("contextmenu",function(d,i){
            d3.event.preventDefault();
            var currColor = this.getAttribute("fill");
            if (currColor != selectedTextColor){
                this.setAttribute("fill",selectedTextColor);
            }else{
                this.setAttribute("fill",null);
            }
        });

      columnText.append("title")
          .text(function(d,i){ return geneNames[i] +": "+ geneDescriptions[i];});
 
  function pair(pair) {
    var cell = d3.select(this).selectAll(".gene")
        .data(pair)
      .enter().append("rect")
        .attr("class", "gene")
        /*
        .attr("x", function(d,i) { return x_axis_scale(i); })
        .attr("width", x_axis_scale.rangeBand())
        .attr("height", x_axis_scale.rangeBand())
        */
        .attr("x", function(d,i) { return i*size; })
        .attr("width", size)
        .attr("height", size)
        .attr("title",function(d){return d[3]})
        .style("fill", 
            function(d) { 
                if(parseFloat(d[3]) < 0){
                    return '#000'
                }else{
                    return color_scale(-Math.log10(parseFloat(d[3]))) 
                }
            })
        .on("click",function(d){
            //d3.json("/detail?GWAS="+d[0]+"&eQTL="+d[1]+"&gene="+d[2],detailcallback);
            d3.select("#detail").remove();
            chart.append("div")
                 .attr("id","detail");
            d3.json("/detail?gene="+d[2],detailcallback);
            //plotDetail(d);
            }
           )
        .append("title")
            .text(function(d){
                return "pval:" + d[3] + "\nqval:" + d[4] + "\nGWAS:" + d[0] + "\neQTL:" + d[1] + "\ngene:" + d[2];
                
            });
     
  }

  function plotDetail(name){
    var a = 1;
  }
    console.log('you called callback! you know how to get data!');
  };

  var draw_pair_json_text = document.getElementById("draw_pair_json_obj").value;
  var draw_pair_json_obj = JSON.parse(draw_pair_json_text);
  //var filtered_gene_names_json_text = document.getElementById("filtered_gene_names_json_obj").value;
  //var filtered_gene_names_json_obj = JSON.parse(filtered_gene_names_json_text);
  //draw_pair(gene_p_qs_json_obj,filtered_gene_names_json_obj);
  draw_pair(draw_pair_json_obj);
    // Code goes here
}

function addDisease(e){
    selected_value = e.value;
    var disease_input_text = document.getElementById("disease_input_text");
    curr_text = disease_input_text.value

    if (curr_text.trim() == ""){
        disease_input_text.value = selected_value;    
    }
    
    //else if( curr_text.indexOf(selected_value) < 0){  
    // allow duplicate because new value can be part of existing value
    // for example,  curr_text = "EXTREME_HEIGHT" and selected_value = "HEIGHT" 
    else{ 
        disease_input_text.value =  disease_input_text.value + ",    " + selected_value;
    }
    //} 

    /*    
    disease_array = curr_text.trim().split("\\s*,\\s");
    new_disease = selected_value.trim();
    
    if (disease_array.indexOf(new_disease) < 0){
        disease_array.push(new_disease);
    }
    disease_input_text.value = disease_array.join(); 
    */
}
