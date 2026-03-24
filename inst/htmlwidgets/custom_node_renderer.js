window.customNodeRenderer = function({ ctx, x, y, state: { selected, hover }, style, label }) {

  const size = style?.size || 80;
  const fillColor = style?.color || '#FFFFFF';

  const width = size * 2.7;
  const height = size;

  return {
    drawNode: function() {

      ctx.save();
      ctx.fillStyle = fillColor;
      ctx.strokeStyle = '#000000';
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.rect(x - width/2, y - height/2, width, height);
      ctx.fill();
      ctx.stroke();
      ctx.restore();

      ctx.save();
      ctx.font = '14px Arial';
      ctx.fillStyle = '#000000';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText(label || '', x, y - height/2 - 5);
      ctx.restore();
    },

    drawExternalLabel: function() {},

    nodeDimensions: { width: width, height: height }
  };
};