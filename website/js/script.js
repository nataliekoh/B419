$(document).ready(function() {
   'use strict';

    $.getJSON('data/tb_prev_2010.json', function(data) {
        var minValue = 100;
        var maxValue = 0;
        for (var key in data) {
            if (data[key] > maxValue) {
                maxValue = data[key];
            } else if (data[key] < minValue) {
                minValue = data[key];
            }
        }

        var paletteScale = d3.scale.linear()
            .domain([minValue, maxValue])
            .range(['#fee0d2', '#67000d']);

        var coloredData = {};
        for (var key in data) {
            var prevValue = Math.round(data[key] * 1000) / 1000;
            coloredData[key] = {
                prev: prevValue,
                fillColor: paletteScale(prevValue)
            };
        }

        var map = new Datamap({
            element: document.getElementById('prev_map'),
            fills: {
                defaultFill: '#CCC'
            },
            data: coloredData,
            geographyConfig: {
                highlightFillColor: function(geo) {
                    return geo['fillColor'] || '#74c476';
                },
                highlightBorderColor: '#74c476',
                popupTemplate: function(geo, data) {
                    if (!data) {
                        return ['<div class="hoverinfo">',
                            '<strong>', geo.properties.name, '</strong>',
                            '<br>Count: <strong>', 'N/A', '</strong>',
                            '</div>'].join('');
                    }

                    return ['<div class="hoverinfo">',
                        '<strong>', geo.properties.name, '</strong>',
                        '<br>Prevalence: <strong>', data.prev, '%</strong>',
                        '</div>'].join('');
                }
            }
        });

        map.legend();
    });


});