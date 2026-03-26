function report = validateManifoldBranch(branch, orbit, mu, opts)
%VALIDATEMANIFOLDBRANCH Basic validation metrics for one manifold branch.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts,'checkApproachSamples'), opts.checkApproachSamples = 10; end

    C = zeros(size(branch.t));
    for k = 1:numel(branch.t)
        C(k) = astro.cr3bp.jacobiConstant(branch.x(k,:).', mu);
    end

    xRef0 = orbit.x(1,:).';
    d0 = norm(branch.x(1,1:6).' - xRef0);

    nCheck = min(opts.checkApproachSamples, size(branch.x,1));
    dEnd = norm(branch.x(nCheck,1:6).' - xRef0);

    report = struct();
    report.type = branch.type;
    report.C0 = C(1);
    report.CmaxError = max(abs(C - C(1)));
    report.initialDistanceToReference = d0;
    report.distanceAfterFirstSamples  = dEnd;
    report.appearsToDepart = dEnd > d0;
end