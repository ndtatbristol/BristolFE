function fn_display_result(nodes, elements, display_options, varargin)
%SUMMARY
%   Displays mesh or results from 2D model
%USAGE
%   fn_display_result(nodes, elements, display_options) to display mesh OR
%   fn_display_result(nodes, elements, display_options, u, global_matrix_node, global_matrix_dof) to display displaced mesh OR
%   fn_display_result(nodes, elements, display_options, u, global_matrix_node, global_matrix_dof, element_colour) to display displaced mesh with elements coloured (e.g. to show stress)
%INPUTS
%   nodes - n x 2 matrix of nodal coordinates. The row number is the node
%   number; columns 1 and 2 are the x and y coordinates of the node.
%   elements - m x 3 matrix of element nodes. The row number is the element
%   number; columns 1, 2 and 3 are the node numbers of the 3 nodes for each
%   triangular element
%   display_options - structured variable allowing optional plotting properties to
%   be set. See below for defaults. In particular:
%   default_options.node_sets_to_plot - allows specific nodes to be plotted
%   in a particular color. It is a vector of structured variables with
%   fields nd and col. nd is a vector of node indices and col is the color
%   (e.g. 'r') in which nodes in that set will be plotted.
%   disp - q x 1 vector of nodal displacements
%   element_colour - m x 1 vector of values used to colour elements (e.g.
%   could be set equal to stress component)
%OUTPUTS
%   none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin >= 4
    u = varargin{1};
else
    u = [];
end
if nargin == 5
    element_colour = varargin{2};
else
    element_colour = ones(size(elements));
end

default_options.draw_elements = 1;
default_options.element_edge_color = 'k';
default_options.mesh_edge_color = 'r';
default_options.draw_mesh_edges = 1;
default_options.node_sets_to_plot = [];
default_options.scale_factor = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display_options = fn_set_default_fields(display_options, default_options);


if ~isempty(u)
    if isempty(display_options.scale_factor)
        model_size = max([max(nodes(:,1)) - min(nodes(:,1)), max(nodes(:,2)) - min(nodes(:,2))]);
        display_options.scale_factor = 1 / max(max(abs(u))) * model_size / 10;
    end
    nodes = nodes + u * display_options.scale_factor;
end

if display_options.draw_elements
    px = reshape(nodes(elements, 1), size(elements)) .';
    py = reshape(nodes(elements, 2), size(elements)) .';
    col = element_colour;
    hold on;
    h = patch(px, py, col', 'EdgeColor', display_options.element_edge_color);
end;

if display_options.draw_mesh_edges
    %get sorted matrix of all element edges
    ed = fn_get_edges(elements);
    %find edges that only occur once (i.e. they are the free edges)
    free_ed = fn_find_free_edges(ed);
    %plot them
    hold on;
    plot(reshape(nodes(free_ed, 1), size(free_ed))', reshape(nodes(free_ed, 2), size(free_ed))', display_options.mesh_edge_color);
end;

if ~isempty(display_options.node_sets_to_plot)
    hold on;
    for ii = 1:length(display_options.node_sets_to_plot)
        plot(nodes(display_options.node_sets_to_plot(ii).nd, 1), nodes(display_options.node_sets_to_plot(ii).nd, 2), display_options.node_sets_to_plot(ii).col);
    end;
end;

axis equal;
axis off;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ed = fn_get_edges(nds)
ed = zeros(prod(size(nds)), 2);
kk = 1;
for ii = 1:size(nds, 1)
    tmp = nds(ii, :);
    tmp = tmp(find(tmp));
    tmp = [tmp, tmp(1)];
    tmp_ed = [tmp(1:end-1); tmp(2:end)]';
    ed(kk: kk + size(tmp_ed,1) - 1, :) = tmp_ed;
    kk = kk + size(tmp_ed,1);
end;
ed = ed(1:kk-1 , :);
ed = sortrows(sort(ed, 2));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function free_ed = fn_find_free_edges(ed)
%add dummy edges at start and end so subsequent logic is general
ed = [[0, 0]; ed; [0, 0]];
ind = find((ed(2:end-1, 1) ~= ed(1:end-2, 1) | ed(2:end-1, 2) ~= ed(1:end-2, 2)) & (ed(2:end-1, 1) ~= ed(3:end, 1) | ed(2:end-1, 2) ~= ed(3:end, 2)));
free_ed = ed(ind + 1, :);
return;

function new_struct = fn_set_default_fields(old_struct, default_struct);
%USAGE
%	new_struct = fn_set_default_fields(old_struct, default_struct);
%SUMMARY
%	Use to add default fields and values to a structured variable, such as
%	options for a function.
%AUTHOR
%	Paul Wilcox (Dec 2003)
%INPUTS
%	old_struct - original structured variable
%	default_struct - structured variable containing default fields and
%	values
%OUTPUTS
%	new_struct - updated structured variable. All existing fields and values
%	in old_struct will be preserved, but any fields found in default_struct
%	and their values will be added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_struct = old_struct;
default_fieldnames = fieldnames(default_struct);
for ii=1:length(default_fieldnames)
	if ~isfield(new_struct, default_fieldnames{ii})
		new_struct = setfield(new_struct, default_fieldnames{ii}, getfield(default_struct, default_fieldnames{ii}));
	end;
end;
return;
