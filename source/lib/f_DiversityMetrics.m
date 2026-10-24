function FD = f_DiversityMetrics(Abundance_in, Traits, Metric, Parameters, Dim)
%Calculating functional diversity measures
%'Rao' Rao's functional biodiversity (quadratic entoropy).
% Botta-Duk�t, Zolt�n 2005
%'FD'  sum of all pairwise distances.
arguments
    Abundance_in
    Traits
    Metric
    Parameters = []
    Dim = 1  %dimension along which we find the diversity index
end

if Dim == 2
    Abundance_in = Abundance_in';
end

for ia = 1:size(Abundance_in, 2)
    Abundance = Abundance_in(:, ia);
    switch Metric
        case 'Rao'
            if numel(Abundance) >2
                P = Abundance/sum(Abundance);
                D = squareform(pdist(Traits,'euclidean'));
                FD(ia) = P' * D * P;
            else
                FD(ia) = 0;
            end
        case 'RaoDist'  %distance should be defined as a Parameters
            P = Abundance/sum(Abundance);
            FD(ia) = P' * Parameters * P;
        case 'FD'  %FD index of functional diveristy
            D = squareform(pdist(Traits,'euclidean'));
            FD(ia) = sum(D(:));
        case 'FDDist'  %FD index, distances should be defined as a Parameters
            FD(ia) = sum(Parameters(:));
        case 'ConvHull'
            switch size(Traits, 2)
                case 2
                    DT = delaunayTriangulation(Traits(:, 1), Traits(:, 2));
                    [~, FD(ia)] = convexHull(DT);
                case 3
                    DT = delaunayTriangulation(Traits(:, 1), Traits(:, 2), Traits(:, 3));
                    [~, FD(ia)] = convexHull(DT);
                otherwise
                    [~,FD(ia)] = convhulln(Traits);
            end
        case 'Hull90'
            %sort abuandances from large to small
            [~, ind] = sort(Abundance, 'desc');
            Abundance = Abundance(ind);
            Traits = Traits(ind, :);
            CumAbund = cumsum(Abundance);
            Idx90 = find(CumAbund>0.95*CumAbund(end), 1, 'first');
            if isempty(Idx90)
                Idx90 = length(CumAbund);
            else
                if Idx90 < 10
                    Idx90 = 10;
                end
                if Idx90 > length(CumAbund)
                    Idx90 = length(CumAbund);
                end
            end
            %get species which contribute 90% of biomass
            Traits = Traits(1:Idx90, :);
            %find min and max traits for these species
            %get volume
            MinTraits = min(Traits, [], 1);
            MaxTraits = max(Traits, [], 1);
            FD(ia) = mean((MaxTraits - MinTraits));
        case 'Richness'
            FD(ia) = sum(Abundance> 0);
        case 'SimpsonESN'
            SumAbundance = sum(Abundance);
            if SumAbundance > 0
                Abundance = Abundance/SumAbundance;
                FD(ia) = 1/( sum(Abundance .* Abundance) );
            else
                FD(ia) = 0;
            end
    end
end
end
