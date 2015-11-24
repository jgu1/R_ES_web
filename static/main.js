$(function() {
  console.log('jquery is working!Jacot');
  createGraph();
});

function createGraph() {
  var margin = 80;
  var w = 700 - 2 * margin, h = 500 - 2 * margin;
  var svg = d3.select("#chart")
                .append("svg")
                    .attr("width", w + 2 * margin)
                    .attr("height", h + 2 * margin)
                .append("svg:g")
                        .attr("transform", "translate(" + margin + ", " + margin + ")");

  var xscale = d3.scale.linear().range([0, w]);
  var yscale = d3.scale.linear().range([h, 0]);

  var xaxis = d3.svg.axis().scale(xscale).ticks(8);
  var yaxis = d3.svg.axis().scale(yscale).ticks(8).orient("left");

  svg.append("svg:g").attr("class", "x axisyo")
    .attr("transform", "translate(0, " + h + ")");
  svg.append("svg:g").attr("class", "y axis");

  svg.select(".x.axisyo").call(xaxis);
  svg.select(".y.axis").call(yaxis);
  
  var callback = function(data){
    console.log('you called callback! you know how to get data!');
  };

  d3.json("/data", callback);
  

  // Code goes here
}
