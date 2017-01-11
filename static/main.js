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
                            wrapper_div.append("div")
                                .attr("id","Manhattan")

                            var selected_geneNames = columnText._groups[0].map(function(text){
                                var oH = text.outerHTML;
                                var selectedTag = 'fill="' + selectedTextColor+'"'; 
                                if (oH.indexOf(selectedTag)>-1){
                                    var geneName = text.__data__;
                                    return geneName;
                                }
                            });

                            d3.json("/Manhattan?geneNames="+selected_geneNames+"&pairNames="+pairNames,Manhattancallback);
                            }
                        ) ;
    append_Manhattan_checkboxes(wrapper_div);
    var a = 1;  
  }

  function append_Manhattan_checkboxes(parentNode){
    d3.select("#show_unaligned_SNPs").remove();
    d3.select("#txt_show_unaligned_SNPs").remove();
    d3.select("#show_untagged_SNPs") .remove();
    d3.select("#txt_show_untagged_SNPs").remove();

 
    var alignedValue = function(d){return d[4]};
    var taggedValue  = function(d){return d[5]};

    parentNode.append("text")
         .text("unaligned")
         .attr("id","txt_show_unaligned_SNPs");
    parentNode.append("input")
         .attr("type","checkbox")
         .attr("checked",true)
         .attr("id","show_unaligned_SNPs")
         .on("change",function(d){
           // the affected SNPs are all unaligned and hence untagged, 
           // brought them back to display only when showing untagged SNPs
           if (this.checked && document.getElementById("show_untagged_SNPs").checked){  
                d3.selectAll(".recteQTL").filter(function(d) {return !alignedValue(d);}) .style("opacity",0.9);
           }else{
                d3.selectAll(".recteQTL").filter(function(d) {return !alignedValue(d);}).style("opacity",0);
                }

         });

    parentNode.append("text")
         .text("untagged")
         .attr("id","txt_show_untagged_SNPs");
    parentNode.append("input")
         .attr("type","checkbox")
         .attr("checked",true)
         .attr("id","show_untagged_SNPs")
         .on("change",function(d){
           if (this.checked){ 
                d3.selectAll(".recteQTL").filter(function(d) {return !taggedValue(d);}).style("opacity",0.9);
                d3.selectAll(".dotGWAS") .filter(function(d) {return !taggedValue(d);}).style("opacity",1);
           
                // there are two types of untagged SNPs: untagged_aligned, untagged_unaligned
                // when showing unaligned_SNPs, both should be brought back to display
                // BUT when "NOT showing unalinged_SNPs", only untagged_aligned should be brought back to display
                if ( ! document.getElementById("show_unaligned_SNPs").checked){
                    d3.selectAll(".recteQTL").filter(function(d) {return !alignedValue(d);}) .style("opacity",0);
                }                

            
            }else{ 
                d3.selectAll(".recteQTL").filter(function(d) {return !taggedValue(d);}).style("opacity", 0);
                d3.selectAll(".dotGWAS") .filter(function(d) {return !taggedValue(d);}).style("opacity",0);
            }
         });

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

 
    d3.select("#Manhattan").html("");
    var Manhattan = d3.select("#Manhattan");

    var margin = {top: 50, right: 20, bottom: 30, left: 40, between: 50};
    var total_width  = 1000;
    var total_height = 180; 
    var Manhattan_width = total_width - margin.left - margin.right; // should not be null, because at least one gene has SNPs
    var Manhattan_height = total_height - margin.top  - margin.bottom;
    var xScaleMax = 3500000000;
    var yScaleMax = 8;
    var yScaleMax_eQTL = 16;
    var recteQTL_size = 5;

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
    for (var i_pair = 0;i_pair < Manhattan_pairNames.length; i_pair ++){
        var curr_Manhattan_pairName = Manhattan_pairNames[i_pair];
        var curr_GWAS_SNPlist = GWAS_SNPlist_dict[curr_Manhattan_pairName];
        //curr_pair_name_x_y = location_pval_chrom_SNPlist_dict[curr_Manhattan_pairName];
        var gene_SNPlist_dict_for_curr_pair = eQTL_gene_SNPlist_dict[curr_Manhattan_pairName];
        //curr_eQTLlist      = eQTL_SNPlist_dict[curr_Manhattan_pairName];
          
        //var zoom = d3.zoom().scaleExtent([1, 8]).on("zoom", zoomed);

        //Zoom in v4
        var gx1 = d3.scaleLinear()
            .range([0,Manhattan_width]);

        var gx2 = d3.scaleLinear()
         
            .range([0,Manhattan_width]);

        gx1.domain([0, xScaleMax]).nice();
         
        gx2.domain([0, xScaleMax]).nice();

        var zoomX = d3.zoom().scaleExtent([0.2, 2000000]).on("zoom", zoomedX);
        //Zoom in v4

        var svg = Manhattan.append("svg")
                .attr("width", Manhattan_width + margin.left + margin.right)
                .attr("height", margin.top + margin.bottom + (Manhattan_height + margin.between) * (NGENES + 1)) //oneGWAS nGenes
                //.attr("height", Manhattan_height + margin.top + margin.bottom)
              .append("g")
                .call(zoomX)
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
                .attr("class","Manhattan_group");
            
        svg.append("rect")
                .attr("width", Manhattan_width)
                .attr("height", Manhattan_height)
                .style("fill", "none")
                .style("pointer-events", "all");
        svg.append("g")
          .attr("class","Manhattan_pairName")
         .append("text")
          .text(curr_Manhattan_pairName);

        function zoomedX(){
            var tr = d3.event.transform;
            gx2 = tr.rescaleX(gx1);
            var x = gx2.domain()[0];
            var y = gx2.domain()[1];
            console.log("gx2: x,y:"+x+ ' '+ y + ' ')
            console.log("###gx1: x,y:"+gx1.domain()[0]+ ' '+ gx1.domain()[1])
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
 

            
        };

        //Zoom in v4
        var xcoord = function(d){
            return gx2(xValue(d));
        };

        //Zoom in v4


        // draw a black line to indicate the location chromosome starts
        
        svg.selectAll(".chrom_starts").data([]).exit().remove();
        svg.selectAll(".gene_starts").data([]).exit().remove();
        for (var i_coordinates = 0; i_coordinates <= NGENES; i_coordinates ++){
            y_shift = (Manhattan_height + margin.between) * i_coordinates
            var chrom_starts =svg.append("g").selectAll(".chrom_starts")
              .data(chrom_starts_data)
              .enter().append("line")
              .attr("class","chrom_starts")
              .attr("x1",function(d){return xScale(d);})
              .attr("x2",function(d){return xScale(d);})
              .attr("y1",y_shift)
              .attr("y2",y_shift + Manhattan_height) 
              .style("opacity", .2)
              .style("stroke","black");
 
        
            var gene_starts =svg.append("g").selectAll(".gene_starts")
              .data(Manhattan_geneNames)
              .enter().append("rect")
              .attr("class","gene_starts")
              .attr("x",function(d){
                var chromStart_chromEnd = gene_location_dict[d]; 
                var chromStart = chromStart_chromEnd[0];
                var x = 0;
                if (chromStart != null){
                    x = xScale(chromStart);
                }
                return x;
              })
              .attr("y",y_shift)
              .attr("width",function(d){
                var chromStart_chromEnd = gene_location_dict[d]; 
                var chromStart = chromStart_chromEnd[0];
                var chromEnd   = chromStart_chromEnd[1];
                var diff       = chromEnd - chromStart;
                var width = 1;
                if (chromStart != null && chromEnd != null){
                    width = xScale(chromEnd) -xScale(chromStart);
                }
                if (width < 1){ width = 1;}
                return width;
              })
              .attr("height",Manhattan_height)
              .attr("fill",function(d){
                return color(d);
              });
        xAxis_shift = y_shift + Manhattan_height;
        // x-axis
        svg.append("g")
          .attr("class", "xAxis")
          .attr("transform", "translate(0," + xAxis_shift + ")")
          //.attr("transform", "translate(0,)")
          .call(xAxis)
        .append("text")
          .attr("class", "label")
          .attr("x", Manhattan_width)
          .attr("y", -6)
          .style("text-anchor", "end")
          .text("Chrom loc");
  
        if (i_coordinates == 0){
            // y-axis
            svg.append("g")
              .attr("class", "yAxis")
              .attr("transform", "translate(0," + y_shift + ")")
              .call(yAxis)
            .append("text")
              .attr("class", "label")
              .attr("y", 6)
              .attr("dy", ".71em")
              .style("text-anchor", "end")
              .text("-log10(p-value)");
        }else{
            // y-axis
            svg.append("g")
              .attr("class", "yAxis_eQTL")
              .call(yAxis_eQTL)
              .attr("transform", "translate(0," + y_shift + ")")
            .append("text")
              .attr("class", "label")
              .attr("y", 6)
              .attr("dy", ".71em")
              .style("text-anchor", "end")
              .text("eQTL: -log10(p-value)");

        }
        }
        svg.append("g").selectAll(".dotGWAS")
          .data(curr_GWAS_SNPlist)
        .enter().append("circle")
          .attr("class", "dotGWAS")
          .attr("r", 3.5)
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
                       //.style("top", (d3.event.pageY - 60) + "px")
                       .style("top", function(d){
                            return (d3.event.pageY - 60) + "px"})
                       .html(d[0] + "(" + d[1] + ")"
                            + " <br/>pval: " + yValue(d) 
                            );
              })
          .on("mouseout", function(d) {
                  tooltip.transition()
                       .duration(500)
                       .style("opacity", 0);
              });

        for (var i_coordinates = 0; i_coordinates < NGENES; i_coordinates ++){
            var y_shift = (Manhattan_height + margin.between) * (i_coordinates + 1); // begin with 1 shift
            var CURR_DRAWING_GENE = Manhattan_geneNames[i_coordinates]; 
            var SNPlist_for_curr_gene_and_curr_pair = gene_SNPlist_dict_for_curr_pair[CURR_DRAWING_GENE];
            /*
            curr_eQTLlist_for_CURR_DRAWING_GENE = new Array();
            for (var i_eSNP = 0; i_eSNP < curr_eQTLlist.length; i_eSNP ++){
                var curr_eSNP_tuple = curr_eQTLlist[i_eSNP]; 
                var gene = cValue(curr_eSNP_tuple);
                if (gene == CURR_DRAWING_GENE){
                    curr_eQTLlist_for_CURR_DRAWING_GENE.push(curr_eSNP_tuple);
                }
            }
            */
             
            //draw eQTLs
            svg.append("g").selectAll(".recteQTL")
              .data(SNPlist_for_curr_gene_and_curr_pair)
            .enter().append("rect")
              .attr("class", "recteQTL")
              .attr("width", recteQTL_size)
              .attr("height",recteQTL_size)
              .attr("x", xMap_eQTL)
              .attr("y", function(d){
                 var shift_within_coordinate = yMap_eQTL(d);
                 var yval = y_shift + shift_within_coordinate;
                 return (yval);
              })
              .attr("aligned",alignedValue)
              .attr("tagged",taggedValue)
              //.style("fill", function(d) {return "black";})
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
                           //.style("top", (d3.event.pageY - 60) + "px")
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
      }
     // draw legend
      var legend_block_size = 10;
      /*
        var legend_wrapper_g = svg.append("g")
           .attr("height",20 * legend_block_size)
           .attr("width", 10 * legend_block_Size);   
      */
    
      var legend = svg.selectAll(".legend")
          .data(color.domain())
        .enter().append("g")
          .attr("class", "legend")
          //.attr("transform", function(d, i) { return "translate(" + Manhattan_width - 10 * legend_bloxk_size + "," + i * legend_block_size + ")"; });
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
          //.attr("x",legend_block_size)
          .attr("y", legend_block_size -1)
          .style("font-size","9px")
          .style("text-anchor", "start")
          .text(function(d) { return d;})

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
                        chart.append("div")
                            .attr("id","Manhattan")

                        var selected_geneNames = columnText._groups[0].map(function(text){
                            var oH = text.outerHTML;
                            var selectedTag = 'fill="' + selectedTextColor+'"'; 
                            if (oH.indexOf(selectedTag)>-1){
                                var geneName = text.__data__;
                                return geneName;
                            }
                        });
                        
                        d3.json("/Manhattan?geneNames="+selected_geneNames+"&pairNames="+pairNames,Manhattancallback);
                        }
                    ) ;
    append_Manhattan_checkboxes(chart);
                                     
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
    }else if( curr_text.indexOf(selected_value) < 0){
        disease_input_text.value =  disease_input_text.value + ",    " + selected_value;
    } 
}
