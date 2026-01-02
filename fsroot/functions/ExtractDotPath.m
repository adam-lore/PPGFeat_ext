function out = ExtractDotPath(data, path)
  for i = 1:numel(path)
    data = [data.(path{i})];
  end

  out = data;
end