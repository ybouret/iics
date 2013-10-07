all:

clean:
	make -C avalanche veryclean
	make -C ff++ clean
	make -C laponite-hs veryclean
	make -C md-bubbles  veryclean
	make -C sludge veryclean
	make -C sludge2 veryclean
	make -C sludge3 veryclean


