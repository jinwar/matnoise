for ip = 1:20
    plot_avg_tomo(ip);
    figure(40)
	title('');
    filename = sprintf('pics/eikonalavgtomo/phasev_%2d',ip);
    filename(find(filename==' '))='0';
    export_fig(filename,'-transparent','-png','-m2')
end
