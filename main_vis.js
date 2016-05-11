/**
 * Created by robertopreste on 29/04/2016.
 */

function createBarplot(sitename) {
    // Create a barplot based on the properties found in the scraped website

    var margin = {top: 30, right: 40, bottom: 100, left: 50},
        width = 800 - margin.left - margin.right,
        height = 500 - margin.top - margin.bottom;

    var x = d3.scale.ordinal()
        .rangeRoundBands([0, width], .1);

    var y = d3.scale.linear()
        .range([height, 0]);

    var xAxis = d3.svg.axis()
        .scale(x)
        .tickFormat(function (d) {
            return d.split("_")[1];
        })
        .orient("bottom");

    var svg = d3.select("#barchart")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .attr("preserveAspectRatio", "xMinYMin meet")
        .append("g")
        .attr("transform", "translate(" + margin.left + ", " + margin.top + ")");

    d3.csv(sitename + "/tagsFound.csv", type, function (error, data) {

        x.domain(data.map(function (d) { return d["Type"]; }));
        y.domain([0, d3.max(data, function (d) { return d["Value"]; }) + 20]);

        var bar = svg.selectAll("g")
            .data(data)
            .enter().append("g")
            .attr("transform", function (d, i) {
                return "translate(" + x(d["Type"]) + ", 0)";
            });

        bar.append("rect")
            .attr("y", function (d) {
                return y(d["Value"]);
            })
            .attr("height", function (d) {
                return height - y(d["Value"]);
            })
            .attr("width", x.rangeBand())
            .attr("fill", function (d) {
                return colorBarplot(d);
            })
            .on("mouseover", function (d) {
                d3.select(this)
                    .attr("fill", "steelblue");
                d3.select(this.parentNode)
                    .selectAll("text")
                    .attr("fill", "#000")
                    .style("font-weight", "bold");
            })
            .on("mouseout", function (d) {
                d3.select(this)
                    .attr("fill", function (d) {
                        return colorBarplot(d);
                    });
                d3.select(this.parentNode)
                    .selectAll("text")
                    .attr("fill", "#666")
                    .style("font-weight", "normal");
            })
            .on("click", function (d) {
                d3.select("#bubbletitle")
                    .selectAll("h2")
                    .remove();
                d3.select("#bubblechart")
                    .selectAll(".bubble")
                    .remove();
                bubbleMaker(sitename, d["Type"]);
                d3.select("#bubbleinfo")
                    .selectAll("table")
                    .remove();
                infoMaker(sitename, d["Type"]);
            });

        bar.append("text")
            .attr("x", x.rangeBand() / 2)
            .attr("y", function (d) {
                return y(d["Value"]) - 20;
            })
            .attr("dy", ".75em")
            .text(function (d) {
                if (d["Value"] != 0) {
                    return (d["Value"]);
                }
            })
            .attr("font-family", "Helvetica Neue, Helvetica, Arial, sans-serif")
            .attr("font-size", "16px")
            .attr("fill", "#666")
            .attr("text-anchor", "middle");

        svg.append("g")
            .attr("transform", "translate(0, " + height + ")")
            .call(xAxis)
            .selectAll("text")
            .style("text-anchor", "end")
            .attr("font-family", "Helvetica Neue, Helvetica, Arial, sans-serif")
            .attr("font-size", "14px")
            .attr("dx", "-.8em")
            .attr("dy", ".15em")
            .attr("transform", "rotate(-50)");

    });
}

function type(d) {
    d["Value"] = +d["Value"];
    return d;
}

function colorBarplot(d) {
    if (d["Type"].charAt(0) == "e") {
        return "#1B62E8";
    } else if (d["Type"].charAt(0) == "o") {
        return "#2AACFF";
    } else if (d["Type"].charAt(0) == "p") {
        return "#81CDFE";
    } else if (d["Type"].charAt(0) == "t") {
        return "#1DFFDE";
    }
}

function bubbleMaker(sitename, csvfile) {
    // Create a bubblechart based on the specific property type selected in the barplot
    
    var diameter = 300,
        color = d3.scale.category20c();

    var bubble = d3.layout.pack()
        .sort(null)
        .size([diameter, diameter])
        .padding(1.5);

    var svg = d3.select("#bubblechart")
        .append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .attr("class", "bubble")
        .attr("id", csvfile)
        .style("padding", "30px")
        .style("margin", "10px")
        .style("background-color", "#fff");

    d3.select("#bubbletitle")
        .append("h2")
        .html(csvfile.split("_")[0].charAt(0).toUpperCase() + csvfile.split("_")[0].slice(1) + " (" + csvfile.split("_")[1] + " type)")
        //.attr("x", diameter / 3)
        //.attr("y", -10)
        //.style("font-weight", "bold")
        //.style("font-size", "16px")
        .style("text-align", "center");

    d3.csv(sitename + "/" + csvfile + ".csv", function (error, data) {
        data = data.map(function (d) {
            d.value = +d["Count"];
            return d;
        });

        var nodes = bubble.nodes({children:data})
            .filter(function (d) {
                return !d.children;
            });

        var bubbles = svg.append("g")
            .attr("transform", "translate(0,0)")
            .selectAll(".bubble")
            .data(nodes)
            .enter();

        bubbles.append("circle")
            .attr("r", function (d) {
                return d.r;
            })
            .attr("cx", function (d) {
                return d.x;
            })
            .attr("cy", function (d) {
                return d.y;
            })
            .style("fill", function (d) {
                return color(d.value);
            });

        bubbles.insert("text")
            .attr("x", function (d) {
                return d.x;
            })
            .attr("y", function (d) {
                return d.y + 5;
            })
            .attr("text-anchor", "middle")
            .text(function (d) {
                return d["Property"];
            })
            .style({"fill": "black",
                    "font-family": "Helvetica Neue, Helvetica, Arial, sans-serif",
                    "font-size": "12px",
                    "font-weight": "bold"});
    })
}

function infoMaker(sitename, csvfile) {
    // Create a table that shows all the details about the properties shown in the bubblechart

    var table = d3.select("#bubbleinfo")
        .append("table");
    var thead = table.append("thead");
    var tbody = table.append("tbody");
    var columns = ["Property", "Type", "Cardinality", "Guideline", "Contr Vocab", "Count"];

    thead.append("tr")
        .selectAll("th")
        .data(columns)
        .enter()
        .append("th")
        .text(function (columns) {return columns;})
        .style({"font-family": "Helvetica Neue, Helvetica, Arial, sans-serif",
                "font-size": "18px"
        });

    d3.csv(sitename + "/" + csvfile + ".csv", function (error, data) {
        data = data.map(function (d) {
            d.value = +d["Count"];
            return d;
        });

        var rows = tbody.selectAll("tr")
            .data(data)
            .enter()
            .append("tr");


        rows.selectAll("td")
            .data(function (row) {
                return columns.map(function (column) {
                    return {column: column, value: row[column], descr: row["Descr"], props: row["Property"]};
                });
            })
            .enter()
            .append("td")
            .style("font-family", "Helvetica Neue, Helvetica, Arial, sans-serif")
            .html(function (d, i) {
                if (i == 0) {
                    return "<a href='https://schema.org/" + d.props + "' target='_blank'>" + d.value + "</a>";
                }
                return d.value;
            })
            .attr("title", function (d) {
                return d.descr;
            });
    });

    return table;

}

function circleProgress(el, sitename) {
    // Create a visualization of compliance for each Bioschemas type
    
    var colors = {
        'pink': '#E1499A',
        'yellow': '#f0ff08',
        'green': '#47e495',
        "event": "#1B62E8",
        "organization": "#2AACFF",
        "person": "#81CDFE",
        "training": "#1DFFDE"
    };
    
    var eventMean = 0;
    var organMean = 0;
    var persMean = 0;
    var trainMean = 0;

    d3.csv(sitename + "/complRate.csv", function (error, data) {
        data = data.map(function (d) {
            d.minScore = +d["MinScore"];
            d.recScore = +d["RecScore"];
            d.optScore = +d["OptScore"];
            d.meanScore = (d.minScore + d.recScore + d.optScore) / 3;

            if (d["Type"].startsWith("event")) {
                eventMean += d.meanScore;
            } else if (d["Type"].startsWith("organization")) {
                organMean += d.meanScore;
            } else if (d["Type"].startsWith("person")) {
                persMean += d.meanScore;
            } else if (d["Type"].startsWith("training")) {
                trainMean += d.meanScore;
            }

            return d;
        });

        eventMean /= 3;
        organMean /= 3;
        persMean /= 3;
        trainMean /= 3;
        
        var endVal, endPar, color;

        if (el == "event") {
            endVal = eventMean;
            endPar = "#event_radial";
            color = colors.event;
        } else if (el == "organization") {
            endVal = organMean;
            endPar = "#organization_radial";
            color = colors.organization;
        } else if (el == "person") {
            endVal = persMean;
            endPar = "#person_radial";
            color = colors.person;
        } else if (el == "training") {
            endVal = trainMean;
            endPar = "#training_radial";
            color = colors.training;
        }

        var radius = 100;
        var border = 10;
        var padding = 30;
        var startPercent = 0;
        var endPercent = endVal / 100;


        var twoPi = Math.PI * 2;
        var formatPercent = d3.format('.0%');
        var boxSize = (radius + padding) * 2;


        var count = Math.abs((endPercent - startPercent) / 0.01);
        var step = endPercent < startPercent ? -0.01 : 0.01;

        var arc = d3.svg.arc()
            .startAngle(0)
            .innerRadius(radius)
            .outerRadius(radius - border);

        var parent = d3.select(endPar);

        parent.append("h2")
            .text(el);

        var svg = parent.append('svg')
            .attr('width', boxSize)
            .attr('height', boxSize);

        var defs = svg.append('defs');

        var filter = defs.append('filter')
            .attr('id', 'blur');

        filter.append('feGaussianBlur')
            .attr('in', 'SourceGraphic')
            .attr('stdDeviation', '7');

        var g = svg.append('g')
            .attr('transform', 'translate(' + boxSize / 2 + ',' + boxSize / 2 + ')');

        var meter = g.append('g')
            .attr('class', 'progress-meter');

        meter.append('path')
            .attr('class', 'background')
            .attr('fill', '#ccc')
            .attr('fill-opacity', 0.5)
            .attr('d', arc.endAngle(twoPi));

        var foreground = meter.append('path')
            .attr('class', 'foreground')
            .attr('fill', color)
            .attr('fill-opacity', 1)
            .attr('stroke', color)
            .attr('stroke-width', 5)
            .attr('stroke-opacity', 1)
            .attr('filter', 'url(#blur)');

        var front = meter.append('path')
            .attr('class', 'foreground')
            .attr('fill', color)
            .attr('fill-opacity', 1);

        var numberText = meter.append('text')
            .attr('fill', '#666')
            .attr('text-anchor', 'middle')
            .attr('dy', '.35em');

        function updateProgress(progress) {
            foreground.attr('d', arc.endAngle(twoPi * progress));
            front.attr('d', arc.endAngle(twoPi * progress));
            numberText.text(formatPercent(progress));
        }

        var progress = startPercent;

        (function loops() {
            updateProgress(progress);

            if (count > 0) {
                count--;
                progress += step;
                setTimeout(loops, 10);
            }
        })();

    });

}

// Dropdown menu for all the scraped websites

var dropbutton = d3.select(".dropdown-content")
    .append("a");

d3.csv("scrapedWebsites.csv", function (error, data) {
    dropbutton.selectAll("a")
        .data(data)
        .enter()
        .append("a")
        .on("click", function (d) {
            d3.select("#barchart")
                .selectAll("svg")
                .remove();
            d3.select("#event_radial")
                .selectAll("svg")
                .remove();
            d3.select("#organization_radial")
                .selectAll("svg")
                .remove();
            d3.select("#person_radial")
                .selectAll("svg")
                .remove();
            d3.select("#training_radial")
                .selectAll("svg")
                .remove();
            var websiteList = d["website"].split("/");
            websiteList.splice(0, 2);
            var websiteName = websiteList.join("_");
            createBarplot(websiteName);
            document.getElementById("barlegendsvg").classList.remove("hide");
            circleProgress("event", websiteName);
            circleProgress("organization", websiteName);
            circleProgress("person", websiteName);
            circleProgress("training", websiteName);
        })
        .html(function (d) {
            return d["name"];
        });

});