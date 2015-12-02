$(function() {
  console.log('jquery is working!Jacot');
  createGraph();
});

function createGraph() {
   var margin = {top: 120, right: 200, bottom: 10, left: 360};
   var canvas2=d3.select("#detail").append("svg")
           .attr("width", 1200)
           .attr("height", 1000)
   var svg2=canvas2
       .append("g")
       .attr("transform", "translate(" + margin.left + "," + margin.top + ")"); 

  //var x_axis_scale = d3.scale.ordinal().rangeBands([0, w])
 
  //var xscale = d3.scale.linear().range([0, w]);
  //var yscale = d3.scale.linear().range([h, 0]);

  var detailcallback = function(data){
    var matrix = [];
    var pairNames = [];
    var SNPNames = [];
    var colSortOrderDesc;
    var size = 15; 
    colSortOrderDesc = false;
    var gene = data.gene;
    var pair_SNP_dict = data.pair_SNP_dict;
    var curr_pair = []
    for (var pairName in pair_SNP_dict) {
        pairNames.push(pairName);
        var SNPs= pair_SNP_dict[pairName];
        var n = SNPs.length;
        curr_pair = []
        SNPs.forEach(function(SNP){
            curr_pair.push(SNP);
            if (SNPNames.length < n){
                SNPNames.push(SNP[0]);
            }
        });  
        matrix.push(curr_pair); 
    } 
    var max_pval = d3.max(matrix, function(pair) {
                   return d3.max(pair.map(function(SNP){return -Math.log10(parseFloat(SNP[2]))})) 
                   }); 
    var min_pval = d3.min(matrix, function(pair) {
                   return d3.min(pair.map(function(SNP){return -Math.log10(parseFloat(SNP[2]))}));
                   }); 

    var color_scale=d3.scale.linear()
                        .domain([min_pval,max_pval])
                        .range(['red','blue']);

    var height_pair = matrix.length * size;
    var width_pair= matrix[0].length * size;

    svg2.html("");// clear canvas2 before drawing a new one
    svg2.append("rect")
      .attr("class", "background")
      .attr("width", width_pair)
      .attr("height", height_pair);

    var pair = svg2.selectAll(".pair2")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair2")
      .attr("transform", function(d, i) { return "translate(0," + i*size + ")"; })
      .each(pair2);

     pair.append("line")
        .attr("x2",width_pair)
        .style("stroke","#fff");

     pair.append("text")
        .attr("x",-6)
        .attr("y",7)
        .attr("dy", ".32em")
        .attr("text-anchor", "end")
        .text(function(d,i) { return pairNames[i]; });

    svg2.selectAll(".column").data([]).exit().remove();
    var column =svg2.selectAll(".column")
              .data(SNPNames)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*size + ")rotate(-90)"; });
               
    column.append("line")
          .style("stroke","#fff");
 

    column.append("line")
        .attr("x1",-height_pair)
        .style("stroke","#fff");


    column.append("text")
      .attr("x", 6)
      .attr("y", 7)
      .attr("width",size)
      .attr("height",size)
      .attr("dy", ".32em")
      .text(function(d, i) { return SNPNames[i]; }); 
 
    function pair2(pair) {
        var cell = d3.select(this).selectAll(".SNP")
            .data(pair)
          .enter().append("rect")
            .attr("class", "SNP")
            .attr("x", function(d,i) { return i*size; })
            .attr("width", size)
            .attr("height", size)
            .style("fill", function(d) { return color_scale(-Math.log10(parseFloat(d[2]))) });
    }
  a = 1;
  }

  var draw_pair = function(data){
    var matrix = [];
    var pairNames = [];
    var geneNames = [];
    var colSortOrderDesc;
    var size = 15; 
    colSortOrderDesc = false;
    for (var pairName in data) {
        if (data.hasOwnProperty(pairName)) {
            pairNames.push(pairName);
            var genes= data[pairName]
            var n = genes.length
            curr_pair = []
            genes.forEach(function(gene){
                curr_pair.push(gene);
                if (geneNames.length < n){
                    geneNames.push(gene[2]);
                }
            });          
            matrix.push(curr_pair)
        }
    }   
    //x_axis_scale.domain(d3.range(matrix[0].length)) // set x_axis_scale's domain to be number of genes in a pair
    var max_pval = d3.max(matrix, function(pair) {
                   return d3.max(pair.map(function(gene){return -Math.log10(parseFloat(gene[3]))})) 
                   }); 
    var min_pval = d3.min(matrix, function(pair) {
                   return d3.min(pair.map(function(gene){return -Math.log10(parseFloat(gene[3]))}));
                   }); 

    var color_scale=d3.scale.linear()
                        .domain([min_pval,max_pval])
                        .range(['red','blue']);

    var height_pair = matrix.length * size;
    var width_pair= matrix[0].length * size;

     var margin = 120;
  
    var w = 900 - 2 * margin, h = 500 - 2 * margin;
    var svg = d3.select("#chart")
              //.style("border-color","#000")
              //.style("border-style","solid")
                  .append("svg")
                    .attr("width", w + 2 * margin)
                    .attr("height",height_pair + 2 * margin)
                  .append("svg:g")
                        .attr("transform", "translate(" + 3*margin + ", " + margin + ")");

    svg.append("rect")
      .attr("class", "background")
      .attr("width", width_pair)
      .attr("height", height_pair);

    var pair = svg.selectAll(".pair")
      .data(matrix)
    .enter().append("g")
      .attr("class", "pair")
      //.attr("transform", function(d, i) { return "translate(0," + x_axis_scale(i) + ")"; })
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
        .text(function(d,i) { return pairNames[i]; });

    svg.selectAll(".column").data([]).exit().remove();
    var column =svg.selectAll(".column")
              .data(geneNames)
              .enter().append("g")
              .attr("class","column")
              .attr("transform", function(d, i) { return "translate(" + i*size + ")rotate(-90)"; });
               
    column.append("line")
          .style("stroke","#fff");
    
    function getColumnI(i){
      var index=0;
      var log2r=[];
      var cd;
      for(index=0; index<matrix.length; index++){
           genes=matrix[index];
           if(genes.length>i){
                log2r.push(parseFloat(genes[i][3]));
           }else{
                log2r.push(0);
           }
      }
      return log2r;
    }

    function sortbylabel(i,colSortOrderDesc){
       var t = svg.transition().duration(3000);
       var log2r=getColumnI(i);
       var sorted; // sorted is zero-based index
       sorted=d3.range(matrix.length).sort(function(a,b){
        if(colSortOrderDesc == true){
            return log2r[a] - log2r[b];
        }else{
            return log2r[b] - log2r[a];
        }
       });
       t.selectAll(".pair")
        .attr("transform", function(d,i) { return "translate(0,"+sorted.indexOf(i) * size+")";});
  }  


    column.append("line")
        .attr("x1",-w)
        .style("stroke","#fff");


    column.append("text")
      .attr("x", 6)
/*    
      .attr("y", x_axis_scale.rangeBand()/2)
      .attr("width",x_axis_scale.rangeBand())
      .attr("height",x_axis_scale.rangeBand())
*/
      .attr("y", 7)
      .attr("width",size)
      .attr("height",size)
      .attr("dy", ".32em")
      .text(function(d, i) { return geneNames[i]; }) 
      .on("click", function(d,i) { colSortOrderDesc=!colSortOrderDesc;sortbylabel(i,colSortOrderDesc);}) ;
 
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
        .style("fill", function(d) { return color_scale(-Math.log10(parseFloat(d[3]))) })
        .on("click",function(d){
            //d3.json("/detail?GWAS="+d[0]+"&eQTL="+d[1]+"&gene="+d[2],detailcallback);
            d3.json("/detail?gene="+d[2],detailcallback);
            //plotDetail(d);
            }
           );  
  }

  function plotDetail(name){
    var a = 1;
  }
    console.log('you called callback! you know how to get data!');
  };

  var gene_p_qs_json_text = document.getElementById("gene_p_qs_json_obj").value;
  var gene_p_qs_json_obj = JSON.parse(gene_p_qs_json_text);
  draw_pair(gene_p_qs_json_obj);
  // Code goes here
}
