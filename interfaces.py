class IParticle:
    def get_origin(self):
        pass
    
    def get_direction(self):
        pass

class IDetector:
    def detects(self, particle: 'IParticle'):
        pass